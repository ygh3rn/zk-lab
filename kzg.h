#pragma once

#include <mcl/bn_c384_256.h>
#include <mcl/bls12_381.hpp>
#include <vector>
#include <chrono>

using namespace mcl;

class KZG {
public:
    struct SetupParams {
        std::vector<G1> g1_powers;  // [g₁, g₁^τ, g₁^τ², ..., g₁^τⁿ]
        std::vector<G2> g2_powers;  // [g₂, g₂^τ] for verification
        size_t max_degree;
    };
    
    struct Commitment {
        G1 commit;
    };
    
    struct Proof {
        G1 witness;
        Fr evaluation;
    };

    // Setup phase - generates structured reference string
    SetupParams Setup(size_t max_degree);
    
    // Commitment phase - commit to polynomial coefficients
    Commitment Commit(const std::vector<Fr>& coefficients, const SetupParams& params);
    
    // Create witness for evaluation at point z
    Proof CreateWitness(const std::vector<Fr>& coefficients, const Fr& z, const SetupParams& params);
    
    // Verify evaluation proof using pairing
    bool VerifyEval(const Commitment& commitment, const Fr& z, const Proof& proof, const SetupParams& params);
    
    // Utility functions
    Fr evaluate_polynomial(const std::vector<Fr>& coefficients, const Fr& x);
    std::vector<Fr> divide_polynomial(const std::vector<Fr>& dividend, const std::vector<Fr>& divisor);
    std::vector<Fr> divide_polynomial_by_linear(const std::vector<Fr>& dividend, const Fr& z);
};