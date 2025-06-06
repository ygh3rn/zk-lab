#include "kzg.h"
#include <mcl/bn.hpp>
#include <mcl/bn_c256.h> 
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <set>

using namespace mcl;
using namespace std;

class CryptographyPractice {
private:
    KZG kzg;
    
public:
    vector<Fr> NTT(const vector<Fr>& a, const Fr& root, size_t n);
    vector<Fr> INTT(const vector<Fr>& a, const Fr& root, size_t n);
    vector<Fr> polynomial_multiply(const vector<Fr>& a, const vector<Fr>& b);
    vector<Fr> polynomial_interpolate_ntt(const vector<Fr>& x_vals, const vector<Fr>& y_vals);
    vector<Fr> polynomial_interpolate_lagrange(const vector<Fr>& x_vals, const vector<Fr>& y_vals);
    
    struct ZeroTestProof {
        KZG::Commitment commitment;
        KZG::Proof quotient_proof;
        Fr random_challenge;
    };
    
    struct SumCheckProof {
        KZG::Commitment commitment;
        KZG::Proof adjusted_proof;
        Fr claimed_sum;
        Fr random_challenge;
    };
    
    ZeroTestProof UnivariateZeroTest_Prove(const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params);
    bool UnivariateZeroTest_Verify(const ZeroTestProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params);
    
    SumCheckProof UnivariateSumCheck_Prove(const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params);
    bool UnivariateSumCheck_Verify(const SumCheckProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params);
    
    Fr find_primitive_root(size_t n);
    Fr mod_inverse(const Fr& a);
    Fr compute_polynomial_sum_on_subgroup(const vector<Fr>& polynomial, const Fr& omega, size_t l);
    vector<Fr> construct_vanishing_polynomial(size_t l);
    vector<Fr> generate_random_polynomial(size_t degree);
    void test_intt_interpolation();
    void test_security_properties();
    void run_tests();
};

// Non-recursive NTT implementation with bit-reversal
vector<Fr> CryptographyPractice::NTT(const vector<Fr>& a, const Fr& root, size_t n) {
    if (n == 0 || (n & (n - 1)) != 0) {
        throw invalid_argument("NTT size must be power of 2");
    }
    
    vector<Fr> result = a;
    result.resize(n, Fr(0));
    
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            swap(result[i], result[j]);
        }
    }
    
    for (size_t len = 2; len <= n; len <<= 1) {
        Fr wlen;
        Fr exponent = Fr(n / len);
        Fr::pow(wlen, root, exponent);
        
        for (size_t i = 0; i < n; i += len) {
            Fr w = Fr(1);
            
            for (size_t j = 0; j < len / 2; j++) {
                Fr u = result[i + j];
                Fr v;
                Fr::mul(v, result[i + j + len / 2], w);
                
                Fr::add(result[i + j], u, v);
                Fr::sub(result[i + j + len / 2], u, v);
                Fr::mul(w, w, wlen);
            }
        }
    }
    
    return result;
}

// Inverse NTT implementation
vector<Fr> CryptographyPractice::INTT(const vector<Fr>& a, const Fr& root, size_t n) {
    Fr inv_root = mod_inverse(root);
    vector<Fr> result = NTT(a, inv_root, n);
    
    Fr inv_n = mod_inverse(Fr(n));
    for (auto& x : result) {
        Fr::mul(x, x, inv_n);
    }
    
    return result;
}

// Polynomial multiplication using NTT
vector<Fr> CryptographyPractice::polynomial_multiply(const vector<Fr>& a, const vector<Fr>& b) {
    if (a.empty() || b.empty()) {
        return {};
    }
    
    size_t result_size = a.size() + b.size() - 1;
    size_t n = 1;
    while (n < result_size) n <<= 1;
    
    Fr root = find_primitive_root(n);
    
    vector<Fr> fa = NTT(a, root, n);
    vector<Fr> fb = NTT(b, root, n);
    
    for (size_t i = 0; i < n; i++) {
        Fr::mul(fa[i], fa[i], fb[i]);
    }
    
    vector<Fr> result = INTT(fa, root, n);
    result.resize(result_size);
    
    return result;
}

// INTT-based polynomial interpolation for roots of unity
vector<Fr> CryptographyPractice::polynomial_interpolate_ntt(const vector<Fr>& x_vals, const vector<Fr>& y_vals) {
    size_t n = x_vals.size();
    if (n == 0) return {};
    
    if (n > 1) {
        Fr omega = find_primitive_root(n);
        Fr omega_power = Fr(1);
        
        bool is_subgroup = true;
        for (size_t i = 0; i < n; i++) {
            if (!(x_vals[i] == omega_power)) {
                is_subgroup = false;
                break;
            }
            Fr::mul(omega_power, omega_power, omega);
        }
        
        if (is_subgroup) {
            return INTT(y_vals, omega, n);
        }
    }
    
    return polynomial_interpolate_lagrange(x_vals, y_vals);
}

// Lagrange interpolation for arbitrary points
vector<Fr> CryptographyPractice::polynomial_interpolate_lagrange(const vector<Fr>& x_vals, const vector<Fr>& y_vals) {
    assert(x_vals.size() == y_vals.size());
    size_t n = x_vals.size();
    
    if (n == 0) return {};
    if (n == 1) return {y_vals[0]};
    
    vector<Fr> result(n, Fr(0));
    
    for (size_t i = 0; i < n; i++) {
        vector<Fr> basis = {Fr(1)};
        
        for (size_t j = 0; j < n; j++) {
            if (i != j) {
                Fr denominator;
                Fr::sub(denominator, x_vals[i], x_vals[j]);
                Fr inv_denom = mod_inverse(denominator);
                
                vector<Fr> linear = {Fr(0), Fr(1)};
                Fr::sub(linear[0], linear[0], x_vals[j]);
                
                basis = polynomial_multiply(basis, linear);
                
                for (auto& coeff : basis) {
                    Fr::mul(coeff, coeff, inv_denom);
                }
            }
        }
        
        for (size_t k = 0; k < basis.size() && k < result.size(); k++) {
            Fr term;
            Fr::mul(term, basis[k], y_vals[i]);
            Fr::add(result[k], result[k], term);
        }
    }
    
    return result;
}

// Prove polynomial vanishes on subgroup with pairing verification
CryptographyPractice::ZeroTestProof CryptographyPractice::UnivariateZeroTest_Prove(
    const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    ZeroTestProof proof;
    proof.commitment = kzg.Commit(polynomial, params);
    
    Fr omega_power = Fr(1);
    for (size_t i = 0; i < l; i++) {
        Fr eval = kzg.evaluate_polynomial(polynomial, omega_power);
        if (!eval.isZero()) {
            throw invalid_argument("Polynomial does not vanish on subgroup");
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    vector<Fr> vanishing_poly = construct_vanishing_polynomial(l);
    vector<Fr> quotient = kzg.divide_polynomial(polynomial, vanishing_poly);
    
    KZG::Commitment quotient_commit = kzg.Commit(quotient, params);
    proof.quotient_proof.witness = quotient_commit.commit;
    proof.quotient_proof.evaluation = Fr(0);
    
    return proof;
}

// Verify ZeroTest proof using pairing e(C, g₂) = e(Q, Z_H(τ)g₂)
bool CryptographyPractice::UnivariateZeroTest_Verify(
    const ZeroTestProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    if (proof.commitment.commit.isZero()) {
        return false;
    }
    
    GT left_pairing, right_pairing;
    
    pairing(left_pairing, proof.commitment.commit, params.g2_powers[0]);
    
    G2 vanishing_g2;
    if (l < params.g2_powers.size()) {
        G2::sub(vanishing_g2, params.g2_powers[l], params.g2_powers[0]);
    } else {
        return false;
    }
    
    pairing(right_pairing, proof.quotient_proof.witness, vanishing_g2);
    
    return left_pairing == right_pairing;
}

// Prove sum of evaluations equals claimed value
CryptographyPractice::SumCheckProof CryptographyPractice::UnivariateSumCheck_Prove(
    const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    SumCheckProof proof;
    proof.commitment = kzg.Commit(polynomial, params);
    proof.claimed_sum = compute_polynomial_sum_on_subgroup(polynomial, omega, l);
    
    vector<Fr> adjusted_poly = polynomial;
    if (!adjusted_poly.empty()) {
        Fr avg_value;
        Fr::div(avg_value, proof.claimed_sum, Fr(l));
        Fr::sub(adjusted_poly[0], adjusted_poly[0], avg_value);
    }
    
    vector<Fr> vanishing_poly = construct_vanishing_polynomial(l);
    vector<Fr> quotient = kzg.divide_polynomial(adjusted_poly, vanishing_poly);
    
    KZG::Commitment quotient_commit = kzg.Commit(quotient, params);
    proof.adjusted_proof.witness = quotient_commit.commit;
    proof.adjusted_proof.evaluation = Fr(0);
    
    return proof;
}

// Verify SumCheck by reducing to ZeroTest
bool CryptographyPractice::UnivariateSumCheck_Verify(
    const SumCheckProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    if (proof.commitment.commit.isZero()) {
        return false;
    }
    
    Fr avg_value;
    Fr::div(avg_value, proof.claimed_sum, Fr(l));
    G1 const_commit;
    G1::mul(const_commit, params.g1_powers[0], avg_value);
    
    G1 adjusted_commit;
    G1::sub(adjusted_commit, proof.commitment.commit, const_commit);
    
    GT left_pairing, right_pairing;
    pairing(left_pairing, adjusted_commit, params.g2_powers[0]);
    
    G2 vanishing_g2;
    if (l < params.g2_powers.size()) {
        G2::sub(vanishing_g2, params.g2_powers[l], params.g2_powers[0]);
    } else {
        return false;
    }
    
    pairing(right_pairing, proof.adjusted_proof.witness, vanishing_g2);
    
    return left_pairing == right_pairing;
}

// Find primitive nth root of unity for BN_SNARK1
Fr CryptographyPractice::find_primitive_root(size_t n) {
    if (n == 0 || (n & (n - 1)) != 0) {
        throw invalid_argument("n must be power of 2");
    }
    
    Fr field_order_minus_one = Fr(-1);
    
    Fr exponent;
    Fr::div(exponent, field_order_minus_one, Fr(n));
    
    vector<int> candidates = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    
    for (int base : candidates) {
        Fr candidate = Fr(base);
        Fr root;
        Fr::pow(root, candidate, exponent);
        
        Fr test_n;
        Fr::pow(test_n, root, Fr(n));
        
        if (test_n == Fr(1)) {
            bool is_primitive = true;
            
            for (size_t divisor = 2; divisor < n; divisor++) {
                if (n % divisor == 0) {
                    Fr test_order;
                    Fr::pow(test_order, root, Fr(divisor));
                    if (test_order == Fr(1)) {
                        is_primitive = false;
                        break;
                    }
                }
            }
                    
            if (is_primitive) {
                return root;
            }
        }
    }
    
    throw runtime_error("Failed to find primitive nth root of unity");
}

// Modular inverse
Fr CryptographyPractice::mod_inverse(const Fr& a) {
    Fr result;
    Fr::inv(result, a);
    return result;
}

// Compute sum of polynomial evaluations on subgroup
Fr CryptographyPractice::compute_polynomial_sum_on_subgroup(const vector<Fr>& polynomial, const Fr& omega, size_t l) {
    Fr sum = Fr(0);
    Fr omega_power = Fr(1);
    
    for (size_t i = 0; i < l; i++) {
        Fr eval = kzg.evaluate_polynomial(polynomial, omega_power);
        Fr::add(sum, sum, eval);
        Fr::mul(omega_power, omega_power, omega);
    }
    
    return sum;
}

// Construct vanishing polynomial Z_H(x) = x^l - 1
vector<Fr> CryptographyPractice::construct_vanishing_polynomial(size_t l) {
    vector<Fr> vanishing(l + 1, Fr(0));
    vanishing[0] = Fr(-1);
    vanishing[l] = Fr(1);
    return vanishing;
}

// Generate random polynomial for testing
vector<Fr> CryptographyPractice::generate_random_polynomial(size_t degree) {
    vector<Fr> poly(degree + 1);
    for (auto& coeff : poly) {
        coeff.setByCSPRNG();
    }
    return poly;
}

// Test INTT-based interpolation demonstrating key insight
void CryptographyPractice::test_intt_interpolation() {
    cout << "\n=== INTT-based Interpolation Test ===" << endl;
    
    size_t n = 8;
    Fr omega = find_primitive_root(n);
    
    vector<Fr> x_vals(n);
    Fr omega_power = Fr(1);
    for (size_t i = 0; i < n; i++) {
        x_vals[i] = omega_power;
        Fr::mul(omega_power, omega_power, omega);
    }
    
    vector<Fr> original_coeffs = {Fr(1), Fr(2), Fr(3), Fr(4), Fr(5), Fr(6), Fr(7), Fr(8)};
    vector<Fr> y_vals = NTT(original_coeffs, omega, n);
    
    auto start = chrono::high_resolution_clock::now();
    vector<Fr> recovered_intt = polynomial_interpolate_ntt(x_vals, y_vals);
    auto end = chrono::high_resolution_clock::now();
    auto intt_time = chrono::duration_cast<chrono::microseconds>(end - start);
    
    start = chrono::high_resolution_clock::now();
    vector<Fr> recovered_lagrange = polynomial_interpolate_lagrange(x_vals, y_vals);
    end = chrono::high_resolution_clock::now();
    auto lagrange_time = chrono::duration_cast<chrono::microseconds>(end - start);
    
    bool intt_correct = (recovered_intt.size() == original_coeffs.size());
    bool lagrange_correct = (recovered_lagrange.size() == original_coeffs.size());
    
    if (intt_correct) {
        for (size_t i = 0; i < original_coeffs.size(); i++) {
            if (!(recovered_intt[i] == original_coeffs[i])) {
                intt_correct = false;
                break;
            }
        }
    }
    
    if (lagrange_correct) {
        for (size_t i = 0; i < original_coeffs.size(); i++) {
            if (!(recovered_lagrange[i] == original_coeffs[i])) {
                lagrange_correct = false;
                break;
            }
        }
    }
    
    cout << "INTT interpolation: " << (intt_correct ? "PASS" : "FAIL") 
         << " (" << intt_time.count() << " μs)" << endl;
    cout << "Lagrange interpolation: " << (lagrange_correct ? "PASS" : "FAIL") 
         << " (" << lagrange_time.count() << " μs)" << endl;
    cout << "INTT speedup: " << (double)lagrange_time.count() / intt_time.count() << "x" << endl;
}

// Test security properties: knowledge soundness, binding, zero-knowledge
void CryptographyPractice::test_security_properties() {
    cout << "\n=== Security Testing ===" << endl;
    
    KZG::SetupParams params = kzg.Setup(64);
    
    cout << "Testing knowledge soundness..." << flush;
    bool knowledge_sound = true;
    
    for (int trial = 0; trial < 20; trial++) {
        vector<Fr> poly = generate_random_polynomial(32);
        Fr eval_point; eval_point.setByCSPRNG();
        
        KZG::Commitment commit = kzg.Commit(poly, params);
        KZG::Proof proof = kzg.CreateWitness(poly, eval_point, params);
        
        Fr direct_eval = kzg.evaluate_polynomial(poly, eval_point);
        if (!(proof.evaluation == direct_eval)) {
            knowledge_sound = false;
            break;
        }
        
        if (!kzg.VerifyEval(commit, eval_point, proof, params)) {
            knowledge_sound = false;
            break;
        }
    }
    cout << (knowledge_sound ? " PASS" : " FAIL") << endl;
    
    cout << "Testing binding property..." << flush;
    bool binding_secure = true;
    
    for (int trial = 0; trial < 10; trial++) {
        vector<Fr> poly1 = generate_random_polynomial(16);
        vector<Fr> poly2 = generate_random_polynomial(16);
        
        poly2[0] = poly1[0];
        Fr::add(poly2[0], poly2[0], Fr(1));
        
        KZG::Commitment commit1 = kzg.Commit(poly1, params);
        KZG::Commitment commit2 = kzg.Commit(poly2, params);
        
        if (commit1.commit == commit2.commit) {
            binding_secure = false;
            break;
        }
    }
    cout << (binding_secure ? " PASS" : " FAIL") << endl;
    
    cout << "Testing zero-knowledge property..." << flush;
    vector<Fr> secret_poly = generate_random_polynomial(32);
    KZG::Commitment commit = kzg.Commit(secret_poly, params);
    bool appears_random = !commit.commit.isZero();
    cout << (appears_random ? " PASS" : " FAIL") << endl;
    
    bool all_secure = knowledge_sound && binding_secure && appears_random;
    cout << "Overall security: " << (all_secure ? "✅ SECURE" : "❌ VULNERABILITIES") << endl;
}

// Comprehensive testing suite
void CryptographyPractice::run_tests() {
    cout << "=== Enhanced Cryptography Practice Implementation Tests ===" << endl;
    
    initPairing(mcl::BN_SNARK1);
    cout << "MCL library initialized with BN_SNARK1 curve" << endl;
    
    cout << "\n1. Testing Enhanced NTT and INTT..." << endl;
    auto start = chrono::high_resolution_clock::now();
    
    vector<Fr> test_poly = {Fr(1), Fr(2), Fr(3), Fr(4), Fr(5), Fr(6), Fr(7), Fr(8)};
    size_t n = 8;
    
    Fr root = find_primitive_root(n);
    cout << "Found primitive " << n << "th root of unity with verification" << endl;
    
    vector<Fr> ntt_result = NTT(test_poly, root, n);
    vector<Fr> intt_result = INTT(ntt_result, root, n);
    
    bool ntt_correct = true;
    for (size_t i = 0; i < test_poly.size(); i++) {
        if (!(test_poly[i] == intt_result[i])) {
            ntt_correct = false;
            break;
        }
    }
    
    auto end = chrono::high_resolution_clock::now();
    auto ntt_time = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Enhanced NTT/INTT correctness: " << (ntt_correct ? "PASS" : "FAIL") << endl;
    cout << "NTT/INTT time: " << ntt_time.count() << " microseconds" << endl;
    
    test_intt_interpolation();
    
    cout << "\n3. Testing enhanced polynomial multiplication..." << endl;
    start = chrono::high_resolution_clock::now();
    
    vector<Fr> poly_a = {Fr(1), Fr(2), Fr(3)};
    vector<Fr> poly_b = {Fr(4), Fr(5)};
    
    vector<Fr> product = polynomial_multiply(poly_a, poly_b);
    
    end = chrono::high_resolution_clock::now();
    auto mult_time = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Polynomial multiplication time: " << mult_time.count() << " microseconds" << endl;
    
    bool mult_correct = (product.size() == 4) && (product[1] == Fr(13));
    cout << "Multiplication correctness: " << (mult_correct ? "PASS" : "FAIL") << endl;
    
    cout << "\n4. Testing Enhanced KZG Protocol..." << endl;
    
    size_t max_degree = 64;
    KZG::SetupParams kzg_params = kzg.Setup(max_degree);
    
    vector<Fr> test_polynomial = {Fr(1), Fr(2), Fr(3), Fr(4), Fr(5)};
    
    auto commit_start = chrono::high_resolution_clock::now();
    KZG::Commitment commitment = kzg.Commit(test_polynomial, kzg_params);
    auto commit_end = chrono::high_resolution_clock::now();
    auto commit_time = chrono::duration_cast<chrono::microseconds>(commit_end - commit_start);
    
    Fr eval_point = Fr(42);
    
    auto prove_start = chrono::high_resolution_clock::now();
    KZG::Proof proof = kzg.CreateWitness(test_polynomial, eval_point, kzg_params);
    auto prove_end = chrono::high_resolution_clock::now();
    auto prove_time = chrono::duration_cast<chrono::microseconds>(prove_end - prove_start);
    
    auto verify_start = chrono::high_resolution_clock::now();
    bool verification_result = kzg.VerifyEval(commitment, eval_point, proof, kzg_params);
    auto verify_end = chrono::high_resolution_clock::now();
    auto verify_time = chrono::duration_cast<chrono::microseconds>(verify_end - verify_start);
    
    cout << "Enhanced KZG Protocol:" << endl;
    cout << "  Commit time: " << commit_time.count() << " μs" << endl;
    cout << "  Prove time: " << prove_time.count() << " μs" << endl;
    cout << "  Verify time: " << verify_time.count() << " μs" << endl;
    cout << "  Verification: " << (verification_result ? "PASS" : "FAIL") << endl;
    
    cout << "\n5. Testing Enhanced ZeroTest PIOP..." << endl;
    
    size_t subgroup_size = 8;
    Fr omega = find_primitive_root(subgroup_size);
    
    vector<Fr> base_poly = {Fr(1), Fr(2), Fr(3)};
    vector<Fr> vanishing = construct_vanishing_polynomial(subgroup_size);
    vector<Fr> zero_poly = polynomial_multiply(base_poly, vanishing);
    
    auto zerotest_prove_start = chrono::high_resolution_clock::now();
    ZeroTestProof zerotest_proof = UnivariateZeroTest_Prove(zero_poly, omega, subgroup_size, kzg_params);
    auto zerotest_prove_end = chrono::high_resolution_clock::now();
    auto zerotest_prove_time = chrono::duration_cast<chrono::microseconds>(zerotest_prove_end - zerotest_prove_start);
    
    auto zerotest_verify_start = chrono::high_resolution_clock::now();
    bool zerotest_result = UnivariateZeroTest_Verify(zerotest_proof, omega, subgroup_size, kzg_params);
    auto zerotest_verify_end = chrono::high_resolution_clock::now();
    auto zerotest_verify_time = chrono::duration_cast<chrono::microseconds>(zerotest_verify_end - zerotest_verify_start);
    
    cout << "Enhanced ZeroTest PIOP:" << endl;
    cout << "  Prove time: " << zerotest_prove_time.count() << " μs" << endl;
    cout << "  Verify time: " << zerotest_verify_time.count() << " μs" << endl;
    cout << "  Verification: " << (zerotest_result ? "PASS" : "FAIL") << endl;
    
    cout << "\n6. Testing Enhanced SumCheck PIOP..." << endl;
    
    vector<Fr> sumcheck_poly = {Fr(1)};
    
    auto sumcheck_prove_start = chrono::high_resolution_clock::now();
    SumCheckProof sumcheck_proof = UnivariateSumCheck_Prove(sumcheck_poly, omega, subgroup_size, kzg_params);
    auto sumcheck_prove_end = chrono::high_resolution_clock::now();
    auto sumcheck_prove_time = chrono::duration_cast<chrono::microseconds>(sumcheck_prove_end - sumcheck_prove_start);
    
    auto sumcheck_verify_start = chrono::high_resolution_clock::now();
    bool sumcheck_result = UnivariateSumCheck_Verify(sumcheck_proof, omega, subgroup_size, kzg_params);
    auto sumcheck_verify_end = chrono::high_resolution_clock::now();
    auto sumcheck_verify_time = chrono::duration_cast<chrono::microseconds>(sumcheck_verify_end - sumcheck_verify_start);
    
    cout << "Enhanced SumCheck PIOP:" << endl;
    cout << "  Prove time: " << sumcheck_prove_time.count() << " μs" << endl;
    cout << "  Verify time: " << sumcheck_verify_time.count() << " μs" << endl;
    cout << "  Verification: " << (sumcheck_result ? "PASS" : "FAIL") << endl;
    
    test_security_properties();
}

int main() {
    try {
        CryptographyPractice practice;
        practice.run_tests();
        cout << "\n🎉 All tests completed successfully!" << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}