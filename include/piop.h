#ifndef PIOP_H
#define PIOP_H

#include <vector>
#include <mcl/bn.hpp>
#include "polynomial.h"
#include "kzg.h"

using namespace mcl::bn;

class UnivariateZeroTest {
private:
    KZG& kzg;
    std::vector<Fr> subgroup;  // H = {ω⁰, ω¹, ..., ωⁿ⁻¹}
    Fr omega;                  // generator of subgroup H
    size_t n;                  // |H|
    Polynomial vanishing_poly; // Z_H(x) = x^n - 1
    
public:
    UnivariateZeroTest(KZG& kzg_instance, size_t subgroup_size);
    
    // CORRECTED: O(1) proof size as required
    struct Proof {
        G1 polynomial_commitment;    // Commitment to f(x)
        G1 quotient_commitment;      // Commitment to q(x) where f(x) = q(x) * Z_H(x)
        // Removed: std::vector<Fr> evaluations;  // This violated O(1) proof size
        // Removed: std::vector<G1> witnesses;    // This violated O(1) proof size
    };
    
    // Prover: O(D)G + O(D)F complexity
    Proof prove(const Polynomial& poly);
    
    // Verifier: O(1)G + O(1)F complexity  
    bool verify(const Proof& proof);
    
private:
    void setup_subgroup();
    Polynomial compute_vanishing_polynomial();  // Local method
    bool check_polynomial_zero_on_subgroup(const Polynomial& poly);
};

class UnivariateSumCheck {
private:
    KZG& kzg;
    std::vector<Fr> subgroup;  // H = {ω⁰, ω¹, ..., ωⁿ⁻¹}
    Fr omega;                  // generator of subgroup H
    size_t n;                  // |H|
    Polynomial vanishing_poly; // Z_H(x) = x^n - 1
    
public:
    UnivariateSumCheck(KZG& kzg_instance, size_t subgroup_size);
    
    // CORRECTED: O(1) proof size as required
    struct Proof {
        G1 polynomial_commitment;    // Commitment to f(x)
        G1 quotient_commitment;      // Commitment to h*(x) 
        Fr claimed_sum;              // ∑_{a∈H} f(a)
        Fr random_challenge;         // Random challenge r used in verification
        Fr quotient_eval;            // h*(r) for verification
        G1 quotient_witness;         // Witness for h*(r)
    };
    
    // Prover: O(D)G + O(D)F complexity
    Proof prove(const Polynomial& poly);
    
    // Verifier: O(1)G + O(1)F complexity
    bool verify(const Proof& proof, const Fr& expected_sum);
    
private:
    void setup_subgroup();
    Polynomial compute_vanishing_polynomial();  // Local method
    Fr compute_sum_over_subgroup(const Polynomial& poly);
    
    // Helper for the univariate sum-check protocol (Lemma 10.2 from Thaler)
    Polynomial create_sum_check_polynomial(const Polynomial& poly, const Fr& claimed_sum);
};

#endif