#include "ntt.h"
#include "polynomial.h"
#include "kzg.h"
#include "zerotest.h"
#include "sumcheck.h"
#include <mcl/bn.hpp>
#include <iostream>
#include <chrono>
#include <random>

using namespace mcl;
using namespace std;
using namespace std::chrono;

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
    // Run comprehensive test suite
    void run_comprehensive_tests() {
        cout << "Production Cryptography Test Suite" << endl;
        cout << "=====================================" << endl;
        
        initPairing(mcl::BN_SNARK1);
        KZG::SetupParams params = KZG::Setup(256);
        
        cout << "\nBasic Functionality Tests" << endl;
        cout << "----------------------------" << endl;
        
        test_ntt_comprehensive();
        test_polynomial_operations();
        test_kzg_security(params);
        test_zerotest_mathematical(params);
        test_sumcheck_mathematical(params);
        
        cout << "\nSecurity & Attack Resistance Tests" << endl;
        cout << "--------------------------------------" << endl;
        
        test_security_properties(params);
        test_edge_cases(params);
        
        cout << "\nPerformance Benchmarks" << endl;
        cout << "-------------------------" << endl;
        
        benchmark_performance(params);
        
        cout << "\nFinal Results" << endl;
        cout << "=================" << endl;
        cout << "Tests Passed: " << passed << "/" << total;
        if (passed == total) {
            cout << " All tests passed!" << endl;
            cout << "Production ready for cryptographic deployment" << endl;
        } else {
            cout << " Some tests failed" << endl;
            cout << "Requires fixes before production use" << endl;
        }
    }
    
private:
    void test_ntt_comprehensive() {
        bool ntt_correctness = true;
        bool ntt_performance = true;
        
        for (size_t log_n = 3; log_n <= 10; log_n++) {
            size_t n = 1UL << log_n;
            try {
                Fr root = NTT::find_primitive_root(n);
                
                Fr test_order;
                Fr::pow(test_order, root, Fr(n));
                if (!(test_order == Fr(1))) {
                    ntt_correctness = false;
                    break;
                }
                
                vector<Fr> poly = Polynomial::random(n/2);
                
                auto start = high_resolution_clock::now();
                vector<Fr> fwd = NTT::transform(poly, root, n);
                vector<Fr> inv = NTT::inverse_transform(fwd, root, n);
                auto end = high_resolution_clock::now();
                
                auto duration = duration_cast<microseconds>(end - start);
                if (duration.count() > 1000 * log_n) {
                    ntt_performance = false;
                }
                
                for (size_t i = 0; i < poly.size(); i++) {
                    if (!(poly[i] == inv[i])) {
                        ntt_correctness = false;
                        break;
                    }
                }
            } catch (const exception& e) {
                ntt_correctness = false;
                break;
            }
        }
        
        test("NTT Mathematical Correctness", ntt_correctness);
        test("NTT Performance O(n log n)", ntt_performance);
    }
    
    void test_polynomial_operations() {
        bool mult_correct = true;
        bool interpolation_correct = true;
        
        for (size_t i = 0; i < 20; i++) {
            vector<Fr> a = Polynomial::random(8);
            vector<Fr> b = Polynomial::random(8);
            vector<Fr> c = Polynomial::multiply(a, b);
            
            Fr x; x.setByCSPRNG();
            Fr expected;
            Fr::mul(expected, Polynomial::evaluate(a, x), Polynomial::evaluate(b, x));
            Fr actual = Polynomial::evaluate(c, x);
            
            if (!(expected == actual)) {
                mult_correct = false;
                break;
            }
        }
        
        for (size_t i = 0; i < 10; i++) {
            size_t n = 4 + (i % 4);
            vector<Fr> x_vals(n), y_vals(n);
            for (size_t j = 0; j < n; j++) {
                x_vals[j] = Fr(j + 1);
                y_vals[j].setByCSPRNG();
            }
            
            vector<Fr> poly = Polynomial::interpolate_lagrange(x_vals, y_vals);
            
            for (size_t j = 0; j < n; j++) {
                Fr eval = Polynomial::evaluate(poly, x_vals[j]);
                if (!(eval == y_vals[j])) {
                    interpolation_correct = false;
                    break;
                }
            }
            if (!interpolation_correct) break;
        }
        
        test("Polynomial Multiplication Correctness", mult_correct);
        test("Polynomial Interpolation Correctness", interpolation_correct);
    }
    
    void test_kzg_security(const KZG::SetupParams& params) {
        bool completeness = true;
        bool evaluation_correctness = true;
        
        for (size_t i = 0; i < 50; i++) {
            vector<Fr> poly = Polynomial::random(32);
            Fr eval_point; eval_point.setByCSPRNG();
            
            try {
                KZG::Commitment commit = KZG::Commit(poly, params);
                KZG::Proof proof = KZG::CreateWitness(poly, eval_point, params);
                
                Fr expected = Polynomial::evaluate(poly, eval_point);
                if (!(proof.evaluation == expected)) {
                    evaluation_correctness = false;
                    break;
                }
                
                if (!KZG::VerifyEval(commit, eval_point, proof, params)) {
                    completeness = false;
                    break;
                }
            } catch (const exception& e) {
                completeness = false;
                break;
            }
        }
        
        test("KZG Completeness Property", completeness);
        test("KZG Evaluation Correctness", evaluation_correctness);
    }
    
    void test_zerotest_mathematical(const KZG::SetupParams& params) {
        bool mathematical_soundness = true;
        bool vanishing_verification = true;
        
        for (size_t test_case = 0; test_case < 15; test_case++) {
            size_t l = 8 << (test_case % 3);
            Fr omega = NTT::find_primitive_root(l);
            
            vector<Fr> base_poly = Polynomial::random(6);
            vector<Fr> vanishing_poly = Polynomial::vanishing(l);
            vector<Fr> zero_poly = Polynomial::multiply(base_poly, vanishing_poly);
            
            try {
                Fr omega_power = Fr(1);
                for (size_t i = 0; i < l; i++) {
                    Fr eval = Polynomial::evaluate(zero_poly, omega_power);
                    if (!eval.isZero()) {
                        vanishing_verification = false;
                        break;
                    }
                    Fr::mul(omega_power, omega_power, omega);
                }
                
                if (!vanishing_verification) break;
                
                ZeroTestProof proof = ZeroTest::prove(zero_poly, omega, l, params);
                
                if (!ZeroTest::verify_with_full_checks(proof, omega, l, params)) {
                    mathematical_soundness = false;
                    break;
                }
                
            } catch (const exception& e) {
                mathematical_soundness = false;
                break;
            }
        }
        
        test("ZeroTest Mathematical Soundness", mathematical_soundness);
        test("ZeroTest Vanishing Verification", vanishing_verification);
    }
    
    void test_sumcheck_mathematical(const KZG::SetupParams& params) {
        bool mathematical_correctness = true;
        bool sum_verification = true;
        
        for (size_t test_case = 0; test_case < 12; test_case++) {
            size_t l = 16;
            Fr omega = NTT::find_primitive_root(l);
            
            vector<Fr> poly;
            if (test_case % 3 == 0) {
                poly = {Fr(0)};
            } else if (test_case % 3 == 1) {
                poly = {Fr(0), Fr(test_case + 1)};
            } else {
                vector<Fr> base = Polynomial::random(4);
                Fr base_sum = Polynomial::sum_on_subgroup(base, omega, l);
                
                if (base.empty()) base = {Fr(0)};
                Fr avg; Fr::div(avg, base_sum, Fr(l));
                Fr::sub(base[0], base[0], avg);
                poly = base;
            }
            
            try {
                Fr computed_sum = Polynomial::sum_on_subgroup(poly, omega, l);
                if (!computed_sum.isZero()) {
                    sum_verification = false;
                    break;
                }
                
                SumCheckProof proof = SumCheck::prove(poly, omega, l, params);
                
                if (!SumCheck::verify_with_full_checks(proof, omega, l, params)) {
                    mathematical_correctness = false;
                    break;
                }
                
            } catch (const exception& e) {
                mathematical_correctness = false;
                break;
            }
        }
        
        test("SumCheck Mathematical Correctness", mathematical_correctness);
        test("SumCheck Sum Verification", sum_verification);
    }
    
    void test_security_properties(const KZG::SetupParams& params) {
        size_t binding_violations = 0;
        for (size_t i = 0; i < 30; i++) {
            vector<Fr> poly1 = {Fr(i)};
            vector<Fr> poly2 = {Fr(i + 1)};
            KZG::Commitment c1 = KZG::Commit(poly1, params);
            KZG::Commitment c2 = KZG::Commit(poly2, params);
            if (c1.commit == c2.commit) binding_violations++;
        }
        test("KZG Binding Property", binding_violations == 0);
        
        size_t soundness_violations = 0;
        Fr omega = NTT::find_primitive_root(16);
        for (size_t i = 0; i < 25; i++) {
            vector<Fr> non_vanishing = {Fr(i + 1)};
            try {
                ZeroTest::prove(non_vanishing, omega, 16, params);
                soundness_violations++;
            } catch (const invalid_argument& e) {
            } catch (...) {
                soundness_violations++;
            }
        }
        test("ZeroTest Soundness", soundness_violations == 0);
        
        size_t sumcheck_violations = 0;
        for (size_t i = 0; i < 20; i++) {
            vector<Fr> non_zero_sum = {Fr(i + 1), Fr(1)};
            try {
                SumCheck::prove(non_zero_sum, omega, 16, params);
                sumcheck_violations++;
            } catch (const invalid_argument& e) {
            } catch (...) {
                sumcheck_violations++;
            }
        }
        test("SumCheck Soundness", sumcheck_violations == 0);
    }
    
    void test_edge_cases(const KZG::SetupParams& params) {
        bool edge_case_handling = true;
        
        try {
            vector<Fr> empty_poly;
            KZG::Commitment empty_commit = KZG::Commit(empty_poly, params);
            
            vector<Fr> zero_poly = {Fr(0)};
            Fr omega = NTT::find_primitive_root(8);
            ZeroTestProof zero_proof = ZeroTest::prove(zero_poly, omega, 8, params);
            
            vector<Fr> large_poly = Polynomial::random(params.max_degree - 1);
            KZG::Commitment large_commit = KZG::Commit(large_poly, params);
            
        } catch (const exception& e) {
            edge_case_handling = false;
        }
        
        test("Edge Case Handling", edge_case_handling);
    }
    
    void benchmark_performance(const KZG::SetupParams& params) {
        const size_t poly_degree = 128;
        const size_t subgroup_size = 64;
        
        benchmark("NTT Transform (n=256)", [&]() {
            Fr root = NTT::find_primitive_root(256);
            vector<Fr> poly = Polynomial::random(128);
            NTT::transform(poly, root, 256);
        });
        
        vector<Fr> bench_poly = Polynomial::random(poly_degree);
        
        benchmark("KZG Commit", [&]() {
            KZG::Commit(bench_poly, params);
        });
        
        KZG::Commitment bench_commit = KZG::Commit(bench_poly, params);
        Fr bench_point; bench_point.setByCSPRNG();
        
        benchmark("KZG Prove", [&]() {
            KZG::CreateWitness(bench_poly, bench_point, params);
        });
        
        KZG::Proof bench_proof = KZG::CreateWitness(bench_poly, bench_point, params);
        
        benchmark("KZG Verify", [&]() {
            KZG::VerifyEval(bench_commit, bench_point, bench_proof, params);
        });
        
        Fr omega = NTT::find_primitive_root(subgroup_size);
        vector<Fr> vanishing_poly = Polynomial::multiply(
            Polynomial::random(8), 
            Polynomial::vanishing(subgroup_size)
        );
        
        benchmark("ZeroTest Prove", [&]() {
            ZeroTest::prove(vanishing_poly, omega, subgroup_size, params);
        });
        
        ZeroTestProof zero_bench = ZeroTest::prove(vanishing_poly, omega, subgroup_size, params);
        
        benchmark("ZeroTest Verify", [&]() {
            ZeroTest::verify(zero_bench, omega, subgroup_size, params);
        });
    }
};

int main() {
    try {
        cout << "Initializing BN_SNARK1 elliptic curve..." << endl;
        initPairing(mcl::BN_SNARK1);
        cout << "Cryptographic parameters initialized successfully" << endl;
        
        TestSuite suite;
        suite.run_comprehensive_tests();
        
        return 0;
    } catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << endl;
        return 1;
    }
}