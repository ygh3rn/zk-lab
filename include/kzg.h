#pragma once

#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;

class KZG {
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

    // Generate structured reference string
    static SetupParams Setup(size_t max_degree);
    
    // Commit to polynomial coefficients
    static Commitment Commit(const std::vector<Fr>& coefficients, const SetupParams& params);
    
    // Create evaluation proof
    static Proof CreateWitness(const std::vector<Fr>& coefficients, const Fr& z, const SetupParams& params);
    
    // Verify evaluation proof using pairing
    static bool VerifyEval(const Commitment& commitment, const Fr& z, const Proof& proof, const SetupParams& params);
};