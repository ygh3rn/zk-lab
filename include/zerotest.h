#pragma once

#include "kzg.h"
#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;

struct ZeroTestProof {
    KZG::Commitment commitment;
    KZG::Proof quotient_proof;
};

class ZeroTest {
public:
    // Prove polynomial vanishes on subgroup H_l = {1, ω, ω², ..., ω^(l-1)}
    static ZeroTestProof prove(const std::vector<Fr>& polynomial, const Fr& omega, 
                              size_t l, const KZG::SetupParams& params);
    
    // Verify ZeroTest proof using pairing check
    static bool verify(const ZeroTestProof& proof, const Fr& omega, 
                      size_t l, const KZG::SetupParams& params);
};