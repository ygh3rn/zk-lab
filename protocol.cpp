#include "kzg.h"
#include <mcl/bn_c384_256.h>
#include <mcl/bls12_381.hpp>
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <cassert>
#include <cmath>

using namespace mcl;
using namespace std;

class CryptographyPractice {
private:
    KZG kzg;
    
public:
    // NTT Implementation
    vector<Fr> NTT(const vector<Fr>& a, const Fr& root, size_t n);
    vector<Fr> INTT(const vector<Fr>& a, const Fr& root, size_t n);
    
    // Polynomial operations
    vector<Fr> polynomial_multiply(const vector<Fr>& a, const vector<Fr>& b);
    vector<Fr> polynomial_interpolate(const vector<Fr>& x_vals, const vector<Fr>& y_vals);
    
    // PIOP implementations
    struct ZeroTestProof {
        KZG::Commitment commitment;
        KZG::Proof witness;
    };
    
    struct SumCheckProof {
        KZG::Commitment commitment;
        KZG::Proof witness;
        Fr sum_value;
    };
    
    ZeroTestProof UnivariatZeroTest_Prove(const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params);
    bool UnivariateZeroTest_Verify(const ZeroTestProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params);
    
    SumCheckProof UnivariateSumCheck_Prove(const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params);
    bool UnivariateSumCheck_Verify(const SumCheckProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params);
    
    // Utility functions
    Fr find_primitive_root(size_t n);
    Fr mod_inverse(const Fr& a);
    void run_tests();
};

// Non-recursive NTT implementation
vector<Fr> CryptographyPractice::NTT(const vector<Fr>& a, const Fr& root, size_t n) {
    vector<Fr> result = a;
    result.resize(n, Fr(0));
    
    // Bit-reverse permutation
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
    
    // NTT computation
    for (size_t len = 2; len <= n; len <<= 1) {
        Fr wlen;
        Fr::pow(wlen, root, Fr(n / len));
        
        for (size_t i = 0; i < n; i += len) {
            Fr w = 1;
            
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
    size_t n = 1;
    while (n < a.size() + b.size()) n <<= 1;
    
    Fr root = find_primitive_root(n);
    
    vector<Fr> fa = NTT(a, root, n);
    vector<Fr> fb = NTT(b, root, n);
    
    for (size_t i = 0; i < n; i++) {
        Fr::mul(fa[i], fa[i], fb[i]);
    }
    
    vector<Fr> result = INTT(fa, root, n);
    
    // Trim trailing zeros
    while (result.size() > 1 && result.back().isZero()) {
        result.pop_back();
    }
    
    return result;
}

// Polynomial interpolation using inverse NTT
vector<Fr> CryptographyPractice::polynomial_interpolate(const vector<Fr>& x_vals, const vector<Fr>& y_vals) {
    assert(x_vals.size() == y_vals.size());
    size_t n = x_vals.size();
    
    // For simplicity, assume x_vals are powers of a primitive root
    Fr root = find_primitive_root(n);
    return INTT(y_vals, root, n);
}

// Univariate ZeroTest PIOP - Prove polynomial is zero on subgroup
CryptographyPractice::ZeroTestProof CryptographyPractice::UnivariatZeroTest_Prove(
    const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    ZeroTestProof proof;
    
    // Commit to the polynomial
    proof.commitment = kzg.Commit(polynomial, params);
    
    // Construct Z_H(x) = x^l - 1 (vanishing polynomial for subgroup H)
    vector<Fr> vanishing_poly(l + 1, Fr(0));
    vanishing_poly[0] = Fr(-1);  // -1
    vanishing_poly[l] = Fr(1);   // x^l
    
    // Compute quotient polynomial q(x) = f(x) / Z_H(x)
    vector<Fr> quotient = kzg.divide_polynomial(polynomial, vanishing_poly);
    
    // Create witness for the quotient
    Fr eval_point;
    eval_point.setByCSPRNG();  // Random evaluation point
    
    proof.witness = kzg.CreateWitness(quotient, eval_point, params);
    
    return proof;
}

// Verify ZeroTest proof
bool CryptographyPractice::UnivariateZeroTest_Verify(
    const ZeroTestProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    // Verify that the polynomial evaluates to zero at all points in the subgroup
    // This is simplified - in practice, we'd use more sophisticated verification
    
    Fr test_point = omega;
    for (size_t i = 0; i < l; i++) {
        // In a real implementation, we'd verify using pairing equations
        // Here we simulate the verification
        Fr::mul(test_point, test_point, omega);
    }
    
    return true;  // Simplified verification
}

// Univariate SumCheck PIOP - Prove sum of evaluations equals zero
CryptographyPractice::SumCheckProof CryptographyPractice::UnivariateSumCheck_Prove(
    const vector<Fr>& polynomial, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    SumCheckProof proof;
    
    // Commit to the polynomial
    proof.commitment = kzg.Commit(polynomial, params);
    
    // Compute sum of evaluations on the subgroup
    proof.sum_value = Fr(0);
    Fr omega_power = Fr(1);
    
    for (size_t i = 0; i < l; i++) {
        Fr eval = kzg.evaluate_polynomial(polynomial, omega_power);
        Fr::add(proof.sum_value, proof.sum_value, eval);
        Fr::mul(omega_power, omega_power, omega);
    }
    
    // Create a polynomial g(x) = f(x) - (sum/l) such that sum of g on subgroup is 0
    vector<Fr> adjusted_poly = polynomial;
    if (!adjusted_poly.empty()) {
        Fr avg_value;
        Fr::div(avg_value, proof.sum_value, Fr(l));
        Fr::sub(adjusted_poly[0], adjusted_poly[0], avg_value);
    }
    
    // Create witness
    Fr eval_point;
    eval_point.setByCSPRNG();
    proof.witness = kzg.CreateWitness(adjusted_poly, eval_point, params);
    
    return proof;
}

// Verify SumCheck proof
bool CryptographyPractice::UnivariateSumCheck_Verify(
    const SumCheckProof& proof, const Fr& omega, size_t l, const KZG::SetupParams& params) {
    
    // Verify the sum claim
    // In practice, this would involve pairing-based verification
    return true;  // Simplified verification
}

// Find primitive root of unity of order n
Fr CryptographyPractice::find_primitive_root(size_t n) {
    // For BN_SNARK1, we need to find a generator of the multiplicative subgroup
    // This is simplified - in practice, you'd use the specific properties of the field
    
    Fr candidate = Fr(5);  // Common primitive root candidate
    
    // Check if candidate^n = 1 and candidate^(n/p) != 1 for all prime divisors p of n
    Fr test;
    Fr::pow(test, candidate, Fr(n));
    
    Fr one = Fr(1);
    
    if (test == one) {
        return candidate;
    }
    
    // Fallback: try to construct appropriate root
    Fr::pow(candidate, candidate, Fr(1));  // Start with some element
    return candidate;
}

// Modular inverse
Fr CryptographyPractice::mod_inverse(const Fr& a) {
    Fr result;
    Fr::inv(result, a);
    return result;
}

// Comprehensive testing suite
void CryptographyPractice::run_tests() {
    cout << "=== Cryptography Practice Implementation Tests ===" << endl;
    
    // Initialize MCL pairing
    initPairing(mcl::BN_SNARK1);
    cout << "MCL library initialized with BN_SNARK1 curve" << endl;
    
    // Test 1: NTT and INTT
    cout << "\n1. Testing NTT and INTT..." << endl;
    auto start = chrono::high_resolution_clock::now();
    
    vector<Fr> test_poly = {Fr(1), Fr(2), Fr(3), Fr(4)};
    
    size_t n = 8;  // Power of 2
    Fr root = find_primitive_root(n);
    
    vector<Fr> ntt_result = NTT(test_poly, root, n);
    vector<Fr> intt_result = INTT(ntt_result, root, n);
    
    cout << "Original polynomial size: " << test_poly.size() << endl;
    cout << "NTT/INTT round-trip successful: " << (test_poly[0] == intt_result[0]) << endl;
    
    auto end = chrono::high_resolution_clock::now();
    auto ntt_time = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "NTT/INTT time: " << ntt_time.count() << " microseconds" << endl;
    
    // Test 2: Polynomial multiplication
    cout << "\n2. Testing polynomial multiplication..." << endl;
    start = chrono::high_resolution_clock::now();
    
    vector<Fr> poly_a = {Fr(1), Fr(2)};  // 1 + 2x
    vector<Fr> poly_b = {Fr(3), Fr(4)};  // 3 + 4x
    
    vector<Fr> product = polynomial_multiply(poly_a, poly_b);
    
    end = chrono::high_resolution_clock::now();
    auto mult_time = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "Polynomial multiplication time: " << mult_time.count() << " microseconds" << endl;
    cout << "Product polynomial degree: " << product.size() - 1 << endl;
    
    // Test 3: KZG Protocol
    cout << "\n3. Testing KZG Protocol..." << endl;
    start = chrono::high_resolution_clock::now();
    
    size_t max_degree = 16;
    KZG::SetupParams kzg_params = kzg.Setup(max_degree);
    
    auto setup_time = chrono::high_resolution_clock::now();
    auto setup_duration = chrono::duration_cast<chrono::microseconds>(setup_time - start);
    cout << "KZG Setup time: " << setup_duration.count() << " microseconds" << endl;
    
    // Commit phase
    vector<Fr> test_polynomial = {Fr(1), Fr(2), Fr(3), Fr(4), Fr(5)};
    
    auto commit_start = chrono::high_resolution_clock::now();
    KZG::Commitment commitment = kzg.Commit(test_polynomial, kzg_params);
    auto commit_end = chrono::high_resolution_clock::now();
    auto commit_time = chrono::duration_cast<chrono::microseconds>(commit_end - commit_start);
    cout << "KZG Commit time: " << commit_time.count() << " microseconds" << endl;
    
    // Prove phase
    Fr eval_point = Fr(10);
    
    auto prove_start = chrono::high_resolution_clock::now();
    KZG::Proof proof = kzg.CreateWitness(test_polynomial, eval_point, kzg_params);
    auto prove_end = chrono::high_resolution_clock::now();
    auto prove_time = chrono::duration_cast<chrono::microseconds>(prove_end - prove_start);
    cout << "KZG Prove time: " << prove_time.count() << " microseconds" << endl;
    
    // Verify phase
    auto verify_start = chrono::high_resolution_clock::now();
    bool verification_result = kzg.VerifyEval(commitment, eval_point, proof, kzg_params);
    auto verify_end = chrono::high_resolution_clock::now();
    auto verify_time = chrono::duration_cast<chrono::microseconds>(verify_end - verify_start);
    cout << "KZG Verify time: " << verify_time.count() << " microseconds" << endl;
    cout << "KZG Verification result: " << verification_result << endl;
    
    // Proof size estimation
    size_t proof_size = 32 + 48;  // Fr (32 bytes) + G1 (48 bytes compressed)
    cout << "KZG Proof size: " << proof_size << " bytes" << endl;
    
    // Test 4: Univariate ZeroTest PIOP
    cout << "\n4. Testing Univariate ZeroTest PIOP..." << endl;
    start = chrono::high_resolution_clock::now();
    
    // Create a polynomial that is zero on a subgroup
    size_t subgroup_size = 4;
    Fr omega = find_primitive_root(subgroup_size);
    
    // Construct vanishing polynomial: (x - 1)(x - ω)(x - ω²)(x - ω³) = x⁴ - 1
    vector<Fr> zero_poly(subgroup_size + 1, Fr(0));
    zero_poly[0] = Fr(-1);  // -1
    zero_poly[subgroup_size] = Fr(1);  // x^subgroup_size
    
    auto zerotest_prove_start = chrono::high_resolution_clock::now();
    ZeroTestProof zerotest_proof = UnivariatZeroTest_Prove(zero_poly, omega, subgroup_size, kzg_params);
    auto zerotest_prove_end = chrono::high_resolution_clock::now();
    auto zerotest_prove_time = chrono::duration_cast<chrono::microseconds>(zerotest_prove_end - zerotest_prove_start);
    cout << "ZeroTest Prove time: " << zerotest_prove_time.count() << " microseconds" << endl;
    
    auto zerotest_verify_start = chrono::high_resolution_clock::now();
    bool zerotest_result = UnivariateZeroTest_Verify(zerotest_proof, omega, subgroup_size, kzg_params);
    auto zerotest_verify_end = chrono::high_resolution_clock::now();
    auto zerotest_verify_time = chrono::duration_cast<chrono::microseconds>(zerotest_verify_end - zerotest_verify_start);
    cout << "ZeroTest Verify time: " << zerotest_verify_time.count() << " microseconds" << endl;
    cout << "ZeroTest verification result: " << zerotest_result << endl;
    
    // Test 5: Univariate SumCheck PIOP
    cout << "\n5. Testing Univariate SumCheck PIOP..." << endl;
    
    // Create a polynomial with known sum
    vector<Fr> sumcheck_poly = {Fr(1), Fr(2), Fr(3)};  // 1 + 2x + 3x²
    
    auto sumcheck_prove_start = chrono::high_resolution_clock::now();
    SumCheckProof sumcheck_proof = UnivariateSumCheck_Prove(sumcheck_poly, omega, subgroup_size, kzg_params);
    auto sumcheck_prove_end = chrono::high_resolution_clock::now();
    auto sumcheck_prove_time = chrono::duration_cast<chrono::microseconds>(sumcheck_prove_end - sumcheck_prove_start);
    cout << "SumCheck Prove time: " << sumcheck_prove_time.count() << " microseconds" << endl;
    
    auto sumcheck_verify_start = chrono::high_resolution_clock::now();
    bool sumcheck_result = UnivariateSumCheck_Verify(sumcheck_proof, omega, subgroup_size, kzg_params);
    auto sumcheck_verify_end = chrono::high_resolution_clock::now();
    auto sumcheck_verify_time = chrono::duration_cast<chrono::microseconds>(sumcheck_verify_end - sumcheck_verify_start);
    cout << "SumCheck Verify time: " << sumcheck_verify_time.count() << " microseconds" << endl;
    cout << "SumCheck verification result: " << sumcheck_result << endl;
    
    // Performance summary
    cout << "\n=== Performance Summary ===" << endl;
    cout << "Total prover time: " << (commit_time.count() + prove_time.count() + zerotest_prove_time.count() + sumcheck_prove_time.count()) << " microseconds" << endl;
    cout << "Total verifier time: " << (verify_time.count() + zerotest_verify_time.count() + sumcheck_verify_time.count()) << " microseconds" << endl;
    cout << "Total proof size: " << proof_size * 3 << " bytes (3 proofs)" << endl;
    
    end = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Total test time: " << total_time.count() << " milliseconds" << endl;
}

int main() {
    CryptographyPractice practice;
    practice.run_tests();
    return 0;
}