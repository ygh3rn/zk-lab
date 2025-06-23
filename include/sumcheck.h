#pragma once

#include "kzg.h"
#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;
using namespace std;

struct SumCheckProof {
    KZG::Commitment f_commitment;
    KZG::Commitment q_commitment;
    Fr challenge;
    KZG::Proof f_eval_proof;
    KZG::Proof q_eval_proof;
    Fr f_eval;
    Fr q_eval;
    Fr claimed_sum;
};

class SumCheck {
public:
    static pair<KZG::Commitment, KZG::Commitment> commit_phase(
        const vector<Fr>& polynomial, const Fr& omega, size_t l, 
        const KZG::SetupParams& params);
    
    static SumCheckProof prove(const vector<Fr>& polynomial, 
                              const Fr& challenge, const Fr& omega, size_t l,
                              const KZG::SetupParams& params);
    
    static bool verify(const SumCheckProof& proof, const Fr& omega, size_t l,
                      const KZG::SetupParams& params);
    
    static bool verify_with_full_checks(const SumCheckProof& proof, const Fr& omega, 
                                       size_t l, const KZG::SetupParams& params);
};