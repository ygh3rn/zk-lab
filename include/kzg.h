#pragma once

#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;
using namespace std;

class KZG {
public:
    struct SetupParams {
        vector<G1> g1_powers;
        vector<G2> g2_powers;
        size_t max_degree;
    };
    
    struct Commitment {
        G1 commit;
    };
    
    struct Proof {
        G1 witness;
        Fr evaluation;
    };
    
    struct BatchProof {
        G1 witness;
        vector<Fr> evaluation_points;
        vector<Fr> evaluations;
    };

    // Generate structured reference string
    static SetupParams Setup(size_t max_degree);
    
    // Commit to polynomial coefficients
    static Commitment Commit(const vector<Fr>& coefficients, const SetupParams& params);
    
    // Create evaluation proof
    static Proof CreateWitness(const vector<Fr>& coefficients, const Fr& z, const SetupParams& params);
    
    // Verify evaluation proof using pairing
    static bool VerifyEval(const Commitment& commitment, const Fr& z, const Proof& proof, const SetupParams& params);
    
    // Create batch evaluation proof for multiple points
    static BatchProof CreateWitnessBatch(const vector<Fr>& coefficients, 
                                        const vector<Fr>& eval_points, 
                                        const SetupParams& params);
    
    // Verify batch evaluation proof using pairing
    static bool VerifyEvalBatch(const Commitment& commitment, const BatchProof& proof, 
                               const SetupParams& params);
};