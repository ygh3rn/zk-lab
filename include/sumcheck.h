#pragma once

#include "kzg.h"
#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;
using namespace std;

struct SumCheckProof {
    KZG::Commitment commitment;
    KZG::Proof adjusted_proof;
    Fr claimed_sum;
};

class SumCheck {
public:
    // Generate a SumCheck proof
    static SumCheckProof prove(const vector<Fr>& polynomial, const Fr& omega, 
                              size_t l, const KZG::SetupParams& params);
    
    // Verify a SumCheck proof using pairing-based verification
    static bool verify(const SumCheckProof& proof, const Fr& omega, 
                      size_t l, const KZG::SetupParams& params);
    
    // Enhanced verification with additional cryptographic checks
    static bool verify_with_full_checks(const SumCheckProof& proof, const Fr& omega, 
                                       size_t l, const KZG::SetupParams& params);

private:
    // For debugging/testing
    static vector<Fr> compute_quotient_by_x(const vector<Fr>& polynomial);
    static bool verify_primitive_root(const Fr& omega, size_t l);
};