#pragma once

#include <mcl/bn_c384_256.h>
#include <mcl/bls12_381.hpp>
#include <vector>
#include <chrono>

using namespace mcl;

class KZG {
private:
    std::vector<G1> g1_powers;  // [g^1, g^α, g^α^2, ..., g^α^n]
    std::vector<G2> g2_powers;  // [h^1, h^α] for verification
    size_t max_degree;
    Fr secret_alpha;  // Only used during setup, then discarded

public:
    struct SetupParams {
        std::vector<G1> g1_powers;
        std::vector<G2> g2_powers;
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
    
    // Verify evaluation proof
    bool VerifyEval(const Commitment& commitment, const Fr& z, const Proof& proof, const SetupParams& params);
    
    // Utility functions
    Fr evaluate_polynomial(const std::vector<Fr>& coefficients, const Fr& x);
    std::vector<Fr> divide_polynomial(const std::vector<Fr>& dividend, const std::vector<Fr>& divisor);
};