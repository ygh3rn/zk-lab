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
    
    struct Proof {
        G1 polynomial_commitment;    // C = [f(τ)]₁
        G1 quotient_commitment;      // Q = [q(τ)]₁  
        Fr random_point;             // r (chosen by verifier)
        Fr polynomial_eval;          // f(r)
        Fr quotient_eval;            // q(r)
        G1 polynomial_witness;       // witness for f(r)
        G1 quotient_witness;         // witness for q(r)
    };
    
    Proof prove(const Polynomial& poly);
    
    bool verify(const Proof& proof);
    
private:
    void setup_subgroup();
    Polynomial compute_vanishing_polynomial();
    Polynomial compute_quotient_polynomial(const Polynomial& poly);
    bool is_primitive_nth_root(const Fr& candidate, size_t n);
};

class UnivariateSumCheck {
private:
    KZG& kzg;
    std::vector<Fr> subgroup;  // H = {ω⁰, ω¹, ..., ωⁿ⁻¹}
    Fr omega;                  // generator of subgroup H
    size_t n;                  // |H|
    
public:
    UnivariateSumCheck(KZG& kzg_instance, size_t subgroup_size);
    
    struct Proof {
        G1 polynomial_commitment;    // C = [f(τ)]₁
        G1 h_star_commitment;        // H* = [h*(τ)]₁
        G1 f_linear_commitment;      // F = [f_linear(τ)]₁
        Fr claimed_sum;              // claimed sum value
        Fr random_point;             // r (chosen by verifier)
        Fr polynomial_eval;          // f(r)
        Fr h_star_eval;              // h*(r)
        Fr f_linear_eval;            // f_linear(r)
        G1 polynomial_witness;       // witness for f(r)
        G1 h_star_witness;          // witness for h*(r)
        G1 f_linear_witness;        // witness for f_linear(r)
    };
    
    Proof prove(const Polynomial& poly);
    
    bool verify(const Proof& proof, const Fr& expected_sum);
    
private:
    void setup_subgroup();
    Fr compute_sum_over_subgroup(const Polynomial& poly);
    Polynomial compute_vanishing_polynomial();
};

#endif