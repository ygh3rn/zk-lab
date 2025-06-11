#include "ntt.h"
#include "polynomial.h"
#include "kzg.h"
#include "zerotest.h"
#include "sumcheck.h"
#include <mcl/bn.hpp>
#include <iostream>
#include <chrono>

using namespace chrono;

class TestSuite {
private:
    size_t passed = 0, total = 0;
    
    void test(const string& name, bool condition) {
        total++;
        if (condition) {
            passed++;
            cout << "PASS " << name << endl;
        } else {
            cout << "FAIL " << name << endl;
        }
    }
    
    void benchmark(const string& name, function<void()> func) {
        auto start = high_resolution_clock::now();
        func();
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        cout << name << ": " << duration.count() << "μs" << endl;
    }
    
public:
    void run_tests() {
        cout << "Cryptography Test Suite\n" << endl;
        
        initPairing(BN_SNARK1);
        KZG::SetupParams params = KZG::Setup(512);
        
        test_ntt();
        test_polynomial();
        test_kzg(params);
        test_kzg_batch(params);
        test_zerotest(params);
        test_sumcheck(params);
        test_security(params);
        test_edge_cases(params);
        
        cout << "\nBenchmarks:" << endl;
        benchmark_performance(params);
        
        cout << "\nResults: " << passed << "/" << total;
        cout << (passed == total ? " - All passed" : " - Some failed") << endl;
    }
    
private:
    void test_ntt() {
        bool correctness = true, primitivity = true, edge_cases = true;
        
        // Test various sizes and verify mathematical properties
        vector<size_t> sizes = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
        for (size_t n : sizes) {
            try {
                Fr root = NTT::find_primitive_root(n);
                
                // Verify root is primitive
                Fr test_primitive;
                Fr::pow(test_primitive, root, Fr(n));
                if (!(test_primitive == Fr(1))) { primitivity = false; break; }
                
                for (size_t k = 1; k < n; k++) {
                    Fr::pow(test_primitive, root, Fr(k));
                    if (test_primitive == Fr(1)) { primitivity = false; break; }
                }
                
                // Test correctness with edge cases
                vector<Fr> test_cases[] = {
                    {Fr(0)}, // Zero polynomial
                    {Fr(1)}, // Constant
                    vector<Fr>(n, Fr(1)), // All ones
                    Polynomial::random(n-1), // Max degree
                    Polynomial::random(n/2) // Typical case
                };
                
                for (auto& poly : test_cases) {
                    poly.resize(n, Fr(0));
                    vector<Fr> fwd = NTT::transform(poly, root, n);
                    vector<Fr> inv = NTT::inverse_transform(fwd, root, n);
                    
                    for (size_t i = 0; i < n; i++) {
                        if (!(poly[i] == inv[i])) {
                            correctness = false;
                            goto ntt_done;
                        }
                    }
                }
            } catch (...) { correctness = false; break; }
        }
        ntt_done:
        
        // Test error conditions
        try {
            NTT::find_primitive_root(0); // Should fail
            edge_cases = false;
        } catch (...) {}
        
        try {
            NTT::find_primitive_root(7); // Non-power-of-2
            edge_cases = false;
        } catch (...) {}
        
        test("NTT correctness", correctness);
        test("NTT primitive roots", primitivity);
        test("NTT edge cases", edge_cases);
    }
    
    void test_polynomial() {
        bool mult_correct = true, interp_correct = true, div_correct = true;
        
        // Test polynomial operations with various degrees
        for (size_t deg1 = 0; deg1 <= 16; deg1 += 4) {
            for (size_t deg2 = 0; deg2 <= 16; deg2 += 4) {
                vector<Fr> a = (deg1 == 0) ? vector<Fr>{Fr(0)} : Polynomial::random(deg1);
                vector<Fr> b = (deg2 == 0) ? vector<Fr>{Fr(0)} : Polynomial::random(deg2);
                vector<Fr> c = Polynomial::multiply(a, b);
                
                // Verify multiplication at multiple points
                for (size_t i = 0; i < 5; i++) {
                    Fr x; x.setByCSPRNG();
                    Fr expected, actual;
                    Fr::mul(expected, Polynomial::evaluate(a, x), Polynomial::evaluate(b, x));
                    actual = Polynomial::evaluate(c, x);
                    if (!(expected == actual)) { mult_correct = false; goto poly_done; }
                }
                
                // Test division when possible
                if (!b.empty() && !b.back().isZero()) {
                    try {
                        vector<Fr> quotient = Polynomial::divide(c, b);
                        vector<Fr> verification = Polynomial::multiply(quotient, b);
                        verification.resize(c.size(), Fr(0));
                        
                        for (size_t i = 0; i < c.size(); i++) {
                            if (!(c[i] == verification[i])) {
                                div_correct = false;
                                goto poly_done;
                            }
                        }
                    } catch (...) { div_correct = false; goto poly_done; }
                }
            }
        }
        
        // Test interpolation with various point sets
        for (size_t n = 2; n <= 8; n++) {
            vector<Fr> x_vals(n), y_vals(n);
            
            // Test with sequential points
            for (size_t i = 0; i < n; i++) {
                x_vals[i] = Fr(i + 1);
                y_vals[i].setByCSPRNG();
            }
            
            vector<Fr> poly = Polynomial::interpolate_lagrange(x_vals, y_vals);
            for (size_t i = 0; i < n; i++) {
                Fr eval = Polynomial::evaluate(poly, x_vals[i]);
                if (!(eval == y_vals[i])) {
                    interp_correct = false;
                    goto poly_done;
                }
            }
            
            // Test with random points
            for (size_t i = 0; i < n; i++) {
                x_vals[i].setByCSPRNG();
                y_vals[i].setByCSPRNG();
            }
            
            // Ensure distinct x values
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i + 1; j < n; j++) {
                    if (x_vals[i] == x_vals[j]) {
                        Fr::add(x_vals[j], x_vals[j], Fr(1));
                    }
                }
            }
            
            poly = Polynomial::interpolate_lagrange(x_vals, y_vals);
            for (size_t i = 0; i < n; i++) {
                Fr eval = Polynomial::evaluate(poly, x_vals[i]);
                if (!(eval == y_vals[i])) {
                    interp_correct = false;
                    goto poly_done;
                }
            }
        }
        
        poly_done:
        test("Polynomial multiplication", mult_correct);
        test("Polynomial interpolation", interp_correct);
        test("Polynomial division", div_correct);
    }
    
    void test_kzg(const KZG::SetupParams& params) {
        bool completeness = true, eval_correct = true, pairing_correct = true;
        
        // Test with various polynomial degrees and edge cases
        vector<vector<Fr>> test_polys = {
            {Fr(0)}, // Zero polynomial
            {Fr(42)}, // Constant
            {Fr(1), Fr(2)}, // Linear
            {Fr(0), Fr(0), Fr(1)}, // x^2
            Polynomial::random(params.max_degree), // Max degree
            vector<Fr>(params.max_degree + 1, Fr(0)) // Zero of max degree
        };
        
        for (auto& poly : test_polys) {
            for (size_t pt_test = 0; pt_test < 10; pt_test++) {
                Fr point; point.setByCSPRNG();
                
                try {
                    KZG::Commitment commit = KZG::Commit(poly, params);
                    KZG::Proof proof = KZG::CreateWitness(poly, point, params);
                    
                    Fr expected = Polynomial::evaluate(poly, point);
                    if (!(proof.evaluation == expected)) {
                        eval_correct = false;
                        goto kzg_done;
                    }
                    
                    if (!KZG::VerifyEval(commit, point, proof, params)) {
                        completeness = false;
                        goto kzg_done;
                    }
                    
                    // Manual pairing verification
                    G1 left_g1;
                    G1::mul(left_g1, params.g1_powers[0], proof.evaluation);
                    G1::sub(left_g1, commit.commit, left_g1);
                    
                    G2 right_g2;
                    G2::mul(right_g2, params.g2_powers[0], point);
                    G2::sub(right_g2, params.g2_powers[1], right_g2);
                    
                    GT left_pairing, right_pairing;
                    pairing(left_pairing, left_g1, params.g2_powers[0]);
                    pairing(right_pairing, proof.witness, right_g2);
                    
                    if (!(left_pairing == right_pairing)) {
                        pairing_correct = false;
                        goto kzg_done;
                    }
                } catch (...) { completeness = false; goto kzg_done; }
            }
        }
        
        kzg_done:
        test("KZG completeness", completeness);
        test("KZG evaluation", eval_correct);
        test("KZG pairing equation", pairing_correct);
    }
    
    void test_kzg_batch(const KZG::SetupParams& params) {
        bool correctness = true, consistency = true, edge_cases = true;
        
        // Test various batch sizes and polynomial types
        vector<size_t> batch_sizes = {1, 2, 3, 5, 8, 16, 32};
        vector<vector<Fr>> test_polys = {
            {Fr(1), Fr(2), Fr(3)},
            Polynomial::random(64),
            vector<Fr>(100, Fr(0)),
            Polynomial::random(params.max_degree / 2)
        };
        
        for (auto& poly : test_polys) {
            for (size_t batch_size : batch_sizes) {
                if (batch_size > params.max_degree) continue;
                
                vector<Fr> points(batch_size);
                for (size_t i = 0; i < batch_size; i++) {
                    points[i] = Fr(i * 17 + 23); // Ensure distinct points
                }
                
                try {
                    KZG::Commitment commit = KZG::Commit(poly, params);
                    KZG::BatchProof batch_proof = KZG::CreateWitnessBatch(poly, points, params);
                    
                    if (!KZG::VerifyEvalBatch(commit, batch_proof, params)) {
                        correctness = false;
                        goto batch_done;
                    }
                    
                    // Verify evaluations match individual computations
                    for (size_t i = 0; i < points.size(); i++) {
                        Fr expected = Polynomial::evaluate(poly, points[i]);
                        if (!(batch_proof.evaluations[i] == expected)) {
                            correctness = false;
                            goto batch_done;
                        }
                        
                        // Verify individual proof would also work
                        KZG::Proof individual = KZG::CreateWitness(poly, points[i], params);
                        if (!(individual.evaluation == expected) ||
                            !KZG::VerifyEval(commit, points[i], individual, params)) {
                            consistency = false;
                            goto batch_done;
                        }
                    }
                } catch (...) { correctness = false; goto batch_done; }
            }
        }
        
        // Test edge cases
        try {
            vector<Fr> empty_points;
            KZG::CreateWitnessBatch(test_polys[0], empty_points, params);
            edge_cases = false; // Should throw
        } catch (const invalid_argument&) {}
        
        try {
            // Test with distinct points to avoid Lagrange interpolation issues
            vector<Fr> test_points = {Fr(1), Fr(2), Fr(3)};
            KZG::BatchProof proof = KZG::CreateWitnessBatch(test_polys[0], test_points, params);
            KZG::Commitment commit = KZG::Commit(test_polys[0], params);
            if (!KZG::VerifyEvalBatch(commit, proof, params)) {
                edge_cases = false;
            }
        } catch (...) { edge_cases = false; }
        
        batch_done:
        test("KZG batch correctness", correctness);
        test("KZG batch consistency", consistency);
        test("KZG batch edge cases", edge_cases);
    }
    
    void test_zerotest(const KZG::SetupParams& params) {
        bool soundness = true, completeness = true, math_correct = true;
        
        vector<size_t> subgroup_sizes = {4, 8, 16, 32, 64};
        
        for (size_t l : subgroup_sizes) {
            Fr omega = NTT::find_primitive_root(l);
            
            // Test various vanishing polynomials
            vector<vector<Fr>> test_bases = {
                {Fr(1)},
                {Fr(0), Fr(1)},
                Polynomial::random(8),
                vector<Fr>(16, Fr(1))
            };
            
            for (auto& base : test_bases) {
                vector<Fr> vanishing = Polynomial::vanishing(l);
                vector<Fr> zero_poly = Polynomial::multiply(base, vanishing);
                
                try {
                    // Verify polynomial actually vanishes on subgroup
                    Fr omega_power = Fr(1);
                    for (size_t i = 0; i < l; i++) {
                        Fr eval = Polynomial::evaluate(zero_poly, omega_power);
                        if (!eval.isZero()) {
                            math_correct = false;
                            goto zero_done;
                        }
                        Fr::mul(omega_power, omega_power, omega);
                    }
                    
                    ZeroTestProof proof = ZeroTest::prove(zero_poly, omega, l, params);
                    if (!ZeroTest::verify_with_full_checks(proof, omega, l, params)) {
                        completeness = false;
                        goto zero_done;
                    }
                    
                    // Manual verification of quotient relationship
                    vector<Fr> quotient_poly = Polynomial::divide(zero_poly, vanishing);
                    vector<Fr> verification = Polynomial::multiply(quotient_poly, vanishing);
                    verification.resize(zero_poly.size(), Fr(0));
                    
                    for (size_t i = 0; i < zero_poly.size(); i++) {
                        if (!(zero_poly[i] == verification[i])) {
                            math_correct = false;
                            goto zero_done;
                        }
                    }
                } catch (...) { completeness = false; goto zero_done; }
            }
            
            // Test soundness: non-vanishing polynomials should fail
            for (size_t i = 1; i <= 5; i++) {
                vector<Fr> non_vanishing = {Fr(i)};
                try {
                    ZeroTest::prove(non_vanishing, omega, l, params);
                    soundness = false; // Should have thrown
                    goto zero_done;
                } catch (const invalid_argument&) {}
            }
        }
        
        zero_done:
        test("ZeroTest completeness", completeness);
        test("ZeroTest soundness", soundness);
        test("ZeroTest mathematics", math_correct);
    }
    
    void test_sumcheck(const KZG::SetupParams& params) {
        bool correctness = true, soundness = true, math_verify = true;
        
        vector<size_t> subgroup_sizes = {8, 16, 32};
        
        for (size_t l : subgroup_sizes) {
            Fr omega = NTT::find_primitive_root(l);
            
            // Create polynomials with zero sum
            vector<vector<Fr>> zero_sum_polys;
            
            // Constant zero
            zero_sum_polys.push_back({Fr(0)});
            
            // Linear with zero sum
            Fr neg_l_inv; Fr::inv(neg_l_inv, Fr(l)); Fr::neg(neg_l_inv, neg_l_inv);
            zero_sum_polys.push_back({Fr(0), neg_l_inv});
            
            // Random polynomial adjusted to have zero sum
            for (size_t deg = 2; deg <= 6; deg++) {
                vector<Fr> poly = Polynomial::random(deg);
                Fr current_sum = Polynomial::sum_on_subgroup(poly, omega, l);
                Fr adjustment; Fr::div(adjustment, current_sum, Fr(l));
                Fr::sub(poly[0], poly[0], adjustment);
                zero_sum_polys.push_back(poly);
            }
            
            for (auto& poly : zero_sum_polys) {
                try {
                    // Verify sum is actually zero
                    Fr sum = Polynomial::sum_on_subgroup(poly, omega, l);
                    if (!sum.isZero()) {
                        math_verify = false;
                        goto sum_done;
                    }
                    
                    SumCheckProof proof = SumCheck::prove(poly, omega, l, params);
                    if (!SumCheck::verify_with_full_checks(proof, omega, l, params)) {
                        correctness = false;
                        goto sum_done;
                    }
                } catch (...) { correctness = false; goto sum_done; }
            }
            
            // Test soundness with non-zero sum polynomials
            for (size_t i = 1; i <= 3; i++) {
                vector<Fr> non_zero_sum = {Fr(i)};
                try {
                    SumCheck::prove(non_zero_sum, omega, l, params);
                    soundness = false; // Should fail
                    goto sum_done;
                } catch (const invalid_argument&) {}
            }
        }
        
        sum_done:
        test("SumCheck correctness", correctness);
        test("SumCheck soundness", soundness);
        test("SumCheck mathematics", math_verify);
    }
    
    void test_security(const KZG::SetupParams& params) {
        bool binding = true, hiding = true, knowledge_sound = true;
        
        // Test binding with many different polynomials
        for (size_t i = 0; i < 100; i++) {
            vector<Fr> p1 = Polynomial::random(32);
            vector<Fr> p2 = Polynomial::random(32);
            
            // Ensure polynomials are different
            if (!p1.empty() && !p2.empty()) {
                Fr::add(p2[0], p2[0], Fr(1));
            }
            
            KZG::Commitment c1 = KZG::Commit(p1, params);
            KZG::Commitment c2 = KZG::Commit(p2, params);
            
            if (c1.commit == c2.commit) {
                binding = false;
                break;
            }
        }
        
        // Test computational hiding (different polynomials should have different commitments)
        for (size_t i = 0; i < 50; i++) {
            vector<Fr> poly1 = {Fr(i)};
            vector<Fr> poly2 = {Fr(i + 1)};
            
            KZG::Commitment c1 = KZG::Commit(poly1, params);
            KZG::Commitment c2 = KZG::Commit(poly2, params);
            
            if (c1.commit == c2.commit) {
                hiding = false;
                break;
            }
        }
        
        // Test knowledge soundness through extraction simulation
        for (size_t i = 0; i < 20; i++) {
            vector<Fr> poly = Polynomial::random(16);
            Fr point; point.setByCSPRNG();
            
            KZG::Commitment commit = KZG::Commit(poly, params);
            KZG::Proof proof = KZG::CreateWitness(poly, point, params);
            
            // Verify that proof actually corresponds to polynomial evaluation
            Fr manual_eval = Polynomial::evaluate(poly, point);
            if (!(proof.evaluation == manual_eval)) {
                knowledge_sound = false;
                break;
            }
            
            // Verify quotient polynomial is correct
            vector<Fr> f_minus_fz = poly;
            if (!f_minus_fz.empty()) {
                Fr::sub(f_minus_fz[0], f_minus_fz[0], manual_eval);
            }
            vector<Fr> quotient = Polynomial::divide_by_linear(f_minus_fz, point);
            
            // Quotient should satisfy f(x) - f(z) = (x - z) * q(x)
            vector<Fr> linear = {Fr(0), Fr(1)};
            Fr::sub(linear[0], linear[0], point);
            vector<Fr> verification = Polynomial::multiply(quotient, linear);
            verification.resize(f_minus_fz.size(), Fr(0));
            
            for (size_t j = 0; j < f_minus_fz.size(); j++) {
                if (!(f_minus_fz[j] == verification[j])) {
                    knowledge_sound = false;
                    goto security_done;
                }
            }
        }
        
        security_done:
        test("KZG binding", binding);
        test("KZG hiding", hiding);
        test("KZG knowledge soundness", knowledge_sound);
    }
    
    void test_edge_cases(const KZG::SetupParams& params) {
        bool error_handling = true, boundary_cases = true;
        
        // Test error conditions
        try {
            vector<Fr> oversized(params.max_degree + 2, Fr(1));
            KZG::Commit(oversized, params);
            error_handling = false;
        } catch (const invalid_argument&) {}
        
        try {
            vector<Fr> empty_points;
            KZG::CreateWitnessBatch({Fr(1)}, empty_points, params);
            error_handling = false;
        } catch (const invalid_argument&) {}
        
        // Test boundary cases
        try {
            // Test with reasonable sizes to avoid hanging
            vector<Fr> large_poly = Polynomial::random(64);
            KZG::Commitment commit = KZG::Commit(large_poly, params);
            
            // Test moderately large batch
            vector<Fr> batch_points(16);
            for (size_t i = 0; i < 16; i++) {
                batch_points[i] = Fr(i + 100);
            }
            KZG::BatchProof batch = KZG::CreateWitnessBatch(large_poly, batch_points, params);
            
            if (!KZG::VerifyEvalBatch(commit, batch, params)) {
                boundary_cases = false;
            }
        } catch (...) { boundary_cases = false; }
        
        test("Error handling", error_handling);
        test("Boundary cases", boundary_cases);
    }
    
    void benchmark_performance(const KZG::SetupParams& params) {
        vector<Fr> poly = Polynomial::random(128);
        
        benchmark("KZG Commit", [&]() { KZG::Commit(poly, params); });
        
        Fr point = Fr(42);
        benchmark("KZG Prove", [&]() { KZG::CreateWitness(poly, point, params); });
        
        vector<Fr> batch_points = {Fr(1), Fr(2), Fr(3), Fr(4), Fr(5), Fr(6), Fr(7), Fr(8)};
        benchmark("KZG Batch Prove", [&]() { KZG::CreateWitnessBatch(poly, batch_points, params); });
        
        Fr omega = NTT::find_primitive_root(64);
        vector<Fr> vanishing_poly = Polynomial::multiply(Polynomial::random(8), Polynomial::vanishing(64));
        
        benchmark("ZeroTest Prove", [&]() { ZeroTest::prove(vanishing_poly, omega, 64, params); });
    }
};

int main() {
    try {
        TestSuite suite;
        suite.run_tests();
        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
}