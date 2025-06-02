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
    
public:
    UnivariateZeroTest(KZG& kzg_instance, size_t subgroup_size);
    
    // Prover functions
    struct Proof {
        G1 polynomial_commitment;
        G1 quotient_commitment;
        std::vector<Fr> evaluations;
        std::vector<G1> witnesses;
    };
    
    Proof prove(const Polynomial& poly);
    
    // Verifier functions
    bool verify(const Proof& proof);
    
private:
    void setup_subgroup();
    Polynomial compute_vanishing_polynomial();
};

class UnivariateSumCheck {
private:
    KZG& kzg;
    std::vector<Fr> subgroup;  // H = {ω⁰, ω¹, ..., ωⁿ⁻¹}
    Fr omega;                  // generator of subgroup H
    size_t n;                  // |H|
    
public:
    UnivariateSumCheck(KZG& kzg_instance, size_t subgroup_size);
    
    // Prover functions
    struct Proof {
        G1 polynomial_commitment;
        G1 quotient_commitment;
        Fr claimed_sum;
        std::vector<G1> witnesses;
    };
    
    Proof prove(const Polynomial& poly);
    
    // Verifier functions
    bool verify(const Proof& proof, const Fr& expected_sum);
    
private:
    void setup_subgroup();
    Fr compute_sum_over_subgroup(const Polynomial& poly);
};

#endif