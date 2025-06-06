#pragma once

#include <mcl/bn.hpp>
#include <mcl/bn_c256.h> 
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

    // Setup phase - generates structured reference string
    SetupParams Setup(size_t max_degree);
    
    // Commitment phase - commit to polynomial coefficients
    Commitment Commit(const vector<Fr>& coefficients, const SetupParams& params);
    
    // Create witness for evaluation at point z
    Proof CreateWitness(const vector<Fr>& coefficients, const Fr& z, const SetupParams& params);
    
    // Verify evaluation proof using pairing
    bool VerifyEval(const Commitment& commitment, const Fr& z, const Proof& proof, const SetupParams& params);
    
    // Evaluate polynomial at point x using Horner's method
    Fr evaluate_polynomial(const vector<Fr>& coefficients, const Fr& x);
    
    // General polynomial division
    vector<Fr> divide_polynomial(const vector<Fr>& dividend, const vector<Fr>& divisor);
    
    // Optimized division by linear polynomial (x - z)
    vector<Fr> divide_polynomial_by_linear(const vector<Fr>& dividend, const Fr& z);
};