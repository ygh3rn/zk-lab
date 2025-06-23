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
    mt19937 rng{random_device{}()};
    
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
    
    Fr random_challenge() {
        string random_str = to_string(rng());
        Fr result;
        result.setStr(random_str);
        return result;
    }
    
    bool vectors_equal(const vector<Fr>& a, const vector<Fr>& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); i++) {
            if (!(a[i] == b[i])) return false;
        }
        return true;
    }
    
public:
    void run_tests() {
        cout << "Interactive Cryptography Test Suite" << endl;
        cout << "====================================" << endl;
        
        initPairing(BN_SNARK1);
        KZG::SetupParams params = KZG::Setup(512);
        
        cout << "\nCore Functionality Tests" << endl;
        cout << "-------------------------" << endl;
        test_ntt();
        test_polynomial();
        test_kzg(params);
        
        cout << "\nInteractive Protocol Tests" << endl;
        cout << "---------------------------" << endl;
        test_zerotest_interactive(params);
        test_sumcheck_interactive(params);
        
        cout << "\nSecurity & Soundness Tests" << endl;
        cout << "---------------------------" << endl;
        test_protocol_soundness(params);
        test_random_challenge_security(params);
        test_commitment_binding(params);
        
        cout << "\nEdge Cases & Error Handling" << endl;
        cout << "----------------------------" << endl;
        test_edge_cases(params);
        
        cout << "\nPerformance Benchmarks" << endl;
        cout << "-----------------------" << endl;
        benchmark_interactive_protocols(params);
        
        cout << "\nFinal Results" << endl;
        cout << "===============" << endl;
        cout << "Tests Passed: " << passed << "/" << total;
        cout << (passed == total ? " - All passed!" : " - Some failed") << endl;
    }
    
private:
    void test_ntt() {
        bool correctness = true;
        
        vector<size_t> sizes = {8, 16, 32, 64, 128, 256};
        for (size_t n : sizes) {
            try {
                Fr root = NTT::find_primitive_root(n);
                
                Fr test_order;
                Fr::pow(test_order, root, Fr(n));
                if (test_order != Fr(1)) {
                    correctness = false;
                    break;
                }
                
                vector<Fr> input = Polynomial::random(n/2);
                vector<Fr> transformed = NTT::transform(input, root, n);
                vector<Fr> reconstructed = NTT::inverse_transform(transformed, root, n);
                
                bool roundtrip_ok = true;
                for (size_t k = 0; k < input.size(); k++) {
                    if (!(input[k] == reconstructed[k])) {
                        roundtrip_ok = false;
                        break;
                    }
                }
                if (!roundtrip_ok) {
                    correctness = false;
                    break;
                }
            } catch (...) { 
                correctness = false; 
                break; 
            }
        }
        
        test("NTT Implementation", correctness);
    }
    
    void test_polynomial() {
        bool operations = true;
        
        try {
            vector<Fr> a = {Fr(1), Fr(2), Fr(3)};
            vector<Fr> b = {Fr(4), Fr(5)};
            vector<Fr> product = Polynomial::multiply(a, b);

            vector<Fr> expected = {Fr(4), Fr(13), Fr(22), Fr(15)};
            if (!vectors_equal(product, expected)) {
                operations = false;
            }
            
            vector<Fr> dividend = {Fr(1), Fr(0), Fr(0), Fr(0), Fr(-1)};
            vector<Fr> divisor = {Fr(-1), Fr(0), Fr(1)};
            vector<Fr> quotient = Polynomial::divide(dividend, divisor);
            
            vector<Fr> reconstruction = Polynomial::multiply(quotient, divisor);
            if (!vectors_equal(reconstruction, dividend)) {
                operations = false;
            }
        } catch (...) { 
            operations = false; 
        }
        
        test("Polynomial Operations", operations);
    }
    
    void test_kzg(const KZG::SetupParams& params) {
        bool correctness = true;
        
        try {
            vector<Fr> poly = {Fr(1), Fr(2), Fr(3), Fr(4)};
            KZG::Commitment commit = KZG::Commit(poly, params);
            
            Fr point = Fr(42);
            KZG::Proof proof = KZG::CreateWitness(poly, point, params);
            
            bool valid = KZG::VerifyEval(commit, point, proof, params);
            if (!valid) correctness = false;
            
        } catch (...) { 
            correctness = false; 
        }
        
        test("KZG Protocol", correctness);
    }
    
    void test_zerotest_interactive(const KZG::SetupParams& params) {
        bool completeness = true;
        bool mathematical_correctness = true;
        
        for (size_t test_case = 0; test_case < 5; test_case++) { 
            try {
                size_t l = 8;
                Fr omega = NTT::find_primitive_root(l);
                
                vector<Fr> base_poly = {Fr(1), Fr(2)};
                vector<Fr> vanishing_poly_coeff = Polynomial::vanishing(l);
                vector<Fr> vanishing_poly = Polynomial::multiply(base_poly, vanishing_poly_coeff);
                
                auto [f_commit, q_commit] = ZeroTest::commit_phase(vanishing_poly, omega, l, params);
                Fr challenge = random_challenge();
                ZeroTestProof proof = ZeroTest::prove(vanishing_poly, challenge, omega, l, params);
                
                if (!ZeroTest::verify_with_full_checks(proof, omega, l, params)) {
                    completeness = false;
                    break;
                }
                
                if (!(proof.f_commitment.commit == f_commit.commit) ||
                    !(proof.q_commitment.commit == q_commit.commit)) {
                    mathematical_correctness = false;
                    break;
                }
                
            } catch (const exception& e) {
                completeness = false;
                cout << "ZeroTest error: " << e.what() << endl;
                break;
            }
        }
        
        test("ZeroTest Completeness", completeness);
        test("ZeroTest Mathematical Correctness", mathematical_correctness);
    }
    
    void test_sumcheck_interactive(const KZG::SetupParams& params) {
        bool completeness = true;
        bool mathematical_correctness = true;
        
        for (size_t test_case = 0; test_case < 5; test_case++) {
            try {
                size_t l = 8;
                Fr omega = NTT::find_primitive_root(l);
                
                vector<Fr> poly;
                if (test_case == 0) {
                    poly = {Fr(0)};
                } else if (test_case == 1) {
                    poly = {Fr(0), Fr(1)};
                } else if (test_case == 2) {
                    poly = {Fr(0), Fr(0), Fr(1)};
                } else {
                    vector<Fr> inner = {Fr(1), Fr(2)};
                    poly = Polynomial::multiply(inner, {Fr(0), Fr(1)});
                }
                
                auto [f_commit, q_commit] = SumCheck::commit_phase(poly, omega, l, params);
                Fr challenge = random_challenge();
                SumCheckProof proof = SumCheck::prove(poly, challenge, omega, l, params);
                
                if (!SumCheck::verify_with_full_checks(proof, omega, l, params)) {
                    completeness = false;
                    break;
                }
                
                Fr computed_sum = Polynomial::sum_on_subgroup(poly, omega, l);
                if (!computed_sum.isZero()) {
                    mathematical_correctness = false;
                    break;
                }
                
            } catch (const exception& e) {
                completeness = false;
                cout << "SumCheck error: " << e.what() << endl;
                break;
            }
        }
        
        test("SumCheck Completeness", completeness);
        test("SumCheck Mathematical Correctness", mathematical_correctness);
    }
    
    void test_protocol_soundness(const KZG::SetupParams& params) {
        size_t zerotest_violations = 0;
        size_t sumcheck_violations = 0;
        
        Fr omega = NTT::find_primitive_root(16);
        
        for (size_t i = 1; i <= 20; i++) {
            try {
                vector<Fr> non_vanishing = {Fr(i)};
                ZeroTest::commit_phase(non_vanishing, omega, 16, params);
                zerotest_violations++;
            } catch (const invalid_argument&) {
            } catch (...) {
                zerotest_violations++;
            }
        }
        
        for (size_t i = 1; i <= 20; i++) {
            try {
                vector<Fr> non_zero_sum = {Fr(i), Fr(1)};
                SumCheck::commit_phase(non_zero_sum, omega, 16, params);
                sumcheck_violations++;
            } catch (const invalid_argument&) {
            } catch (...) {
                sumcheck_violations++;
            }
        }
        
        test("ZeroTest Soundness", zerotest_violations == 0);
        test("SumCheck Soundness", sumcheck_violations == 0);
    }
    
    void test_random_challenge_security(const KZG::SetupParams& params) {
        bool challenge_independence = true;
        
        try {
            size_t l = 16;
            Fr omega = NTT::find_primitive_root(l);
            
            vector<Fr> base = {Fr(1), Fr(2)};
            vector<Fr> vanishing_coeff = Polynomial::vanishing(l);
            vector<Fr> vanishing_poly = Polynomial::multiply(base, vanishing_coeff);
            
            vector<Fr> challenges;
            vector<ZeroTestProof> proofs;
            
            for (size_t i = 0; i < 5; i++) {
                Fr challenge = random_challenge();
                challenges.push_back(challenge);
                
                ZeroTestProof proof = ZeroTest::prove(vanishing_poly, challenge, omega, l, params);
                proofs.push_back(proof);
                
                if (!ZeroTest::verify_with_full_checks(proof, omega, l, params)) {
                    challenge_independence = false;
                    break;
                }
            }
            
        } catch (...) {
            challenge_independence = false;
        }
        
        test("Random Challenge Security", challenge_independence);
    }
    
    void test_commitment_binding(const KZG::SetupParams& params) {
        bool binding_property = true;
        
        for (size_t i = 0; i < 20; i++) {
            vector<Fr> poly1 = {Fr(i)};
            vector<Fr> poly2 = {Fr(i + 1)};
            
            KZG::Commitment c1 = KZG::Commit(poly1, params);
            KZG::Commitment c2 = KZG::Commit(poly2, params);
            
            if (c1.commit == c2.commit) {
                binding_property = false;
                break;
            }
        }
        
        test("KZG Binding Property", binding_property);
    }
    
    void test_edge_cases(const KZG::SetupParams& params) {
        bool edge_case_handling = true;
        
        try {
            Fr omega = NTT::find_primitive_root(8);
            Fr challenge = random_challenge();
            
            vector<Fr> zero_poly = {Fr(0)};
            try {
                SumCheckProof sum_proof = SumCheck::prove(zero_poly, challenge, omega, 8, params);
                if (!SumCheck::verify_with_full_checks(sum_proof, omega, 8, params)) {
                    edge_case_handling = false;
                }
            } catch (const exception& e) {
                cout << "Zero poly error: " << e.what() << endl;
                edge_case_handling = false;
            }
            
            vector<Fr> x_poly = {Fr(0), Fr(1)};
            try {
                SumCheckProof x_proof = SumCheck::prove(x_poly, challenge, omega, 8, params);
                if (!SumCheck::verify_with_full_checks(x_proof, omega, 8, params)) {
                    edge_case_handling = false;
                }
            } catch (const exception& e) {
                cout << "X poly error: " << e.what() << endl;
                edge_case_handling = false;
            }
            
        } catch (const exception& e) {
            cout << "Edge case error: " << e.what() << endl;
            edge_case_handling = false;
        }
        
        test("Edge Case Handling", edge_case_handling);
    }
    
    void benchmark_interactive_protocols(const KZG::SetupParams& params) {
        size_t l = 64;
        Fr omega = NTT::find_primitive_root(l);
        
        vector<Fr> base_poly = Polynomial::random(8);
        vector<Fr> vanishing_coeff = Polynomial::vanishing(l);
        vector<Fr> vanishing_poly = Polynomial::multiply(base_poly, vanishing_coeff);
        Fr challenge = random_challenge();
        
        benchmark("ZeroTest Commit Phase", [&]() {
            ZeroTest::commit_phase(vanishing_poly, omega, l, params);
        });
        
        benchmark("ZeroTest Prove", [&]() {
            ZeroTest::prove(vanishing_poly, challenge, omega, l, params);
        });
        
        ZeroTestProof zero_proof = ZeroTest::prove(vanishing_poly, challenge, omega, l, params);
        benchmark("ZeroTest Verify", [&]() {
            ZeroTest::verify_with_full_checks(zero_proof, omega, l, params);
        });
        
        vector<Fr> sum_poly = Polynomial::random(8);
        Fr sum = Polynomial::sum_on_subgroup(sum_poly, omega, l);
        Fr adjustment;
        Fr::div(adjustment, sum, Fr(l));
        Fr::sub(sum_poly[0], sum_poly[0], adjustment);
        
        benchmark("SumCheck Commit Phase", [&]() {
            SumCheck::commit_phase(sum_poly, omega, l, params);
        });
        
        benchmark("SumCheck Prove", [&]() {
            SumCheck::prove(sum_poly, challenge, omega, l, params);
        });
        
        SumCheckProof sum_proof = SumCheck::prove(sum_poly, challenge, omega, l, params);
        benchmark("SumCheck Verify", [&]() {
            SumCheck::verify_with_full_checks(sum_proof, omega, l, params);
        });
    }
};

int main() {
    try {
        TestSuite suite;
        suite.run_tests();
        return 0;
    } catch (const exception& e) {
        cerr << "Test suite error: " << e.what() << endl;
        return 1;
    }
}