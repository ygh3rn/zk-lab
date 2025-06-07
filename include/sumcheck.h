#pragma once

#include "kzg.h"
#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;

struct SumCheckProof {
    KZG::Commitment commitment;
    KZG::Proof adjusted_proof;
    Fr claimed_sum;
};

class SumCheck {
public:
    // Prove sum of polynomial evaluations on subgroup H_l equals claimed value
    static SumCheckProof prove(const std::vector<Fr>& polynomial, const Fr& omega, 
                              size_t l, const KZG::SetupParams& params);
    
    // Verify SumCheck proof using pairing check
    static bool verify(const SumCheckProof& proof, const Fr& omega, 
                      size_t l, const KZG::SetupParams& params);
};