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
    
    // Precomputed commitments for pairing-based verification
    G1 vanishing_commitment_g1;  // [Z_H(τ)]_1
    G2 vanishing_commitment_g2;  // [Z_H(τ)]_2
    
public:
    UnivariateZeroTest(KZG& kzg_instance, size_t subgroup_size);
    
    // Prover functions
    struct Proof {
        G1 polynomial_commitment;  // [f(τ)]_1
        G1 quotient_commitment;    // [q(τ)]_1 where f(x) = q(x) · Z_H(x)
    };
    
    Proof prove(const Polynomial& poly);
    
    // Verifier functions - FULL PAIRING-BASED VERIFICATION
    bool verify(const Proof& proof);
    
private:
    void setup_subgroup();
    void precompute_vanishing_commitment();
    Polynomial compute_vanishing_polynomial();
};

class UnivariateSumCheck {
private:
    KZG& kzg;
    std::vector<Fr> subgroup;  // H = {ω⁰, ω¹, ..., ωⁿ⁻¹}
    Fr omega;                  // generator of subgroup H
    size_t n;                  // |H|
    
    // Precomputed commitments for pairing-based verification
    G1 vanishing_commitment_g1;  // [Z_H(τ)]_1
    G2 vanishing_commitment_g2;  // [Z_H(τ)]_2
    
    // Auxiliary polynomial for verification
    Polynomial auxiliary_poly;  // g(x) = f(x) - S/n
    
public:
    UnivariateSumCheck(KZG& kzg_instance, size_t subgroup_size);
    
    // Prover functions
    struct Proof {
        G1 polynomial_commitment;  // [f(τ)]_1
        G1 quotient_commitment;    // [q(τ)]_1 for sum constraint
        Fr claimed_sum;            // S = ∑_{h∈H} f(h)
    };
    
    Proof prove(const Polynomial& poly);
    
    // Verifier functions - FULL PAIRING-BASED VERIFICATION
    bool verify(const Proof& proof, const Fr& expected_sum);
    
private:
    void setup_subgroup();
    void precompute_vanishing_commitment();
    Fr compute_sum_over_subgroup(const Polynomial& poly);
};

#endif