#pragma once

#include "kzg.h"
#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;
using namespace std;

struct ZeroTestProof {
    KZG::Commitment commitment;
    KZG::Proof quotient_proof;
};

class ZeroTest {
public:
    // Generate a ZeroTest proof
    static ZeroTestProof prove(const vector<Fr>& polynomial, const Fr& omega, 
                              size_t l, const KZG::SetupParams& params);
    
    // Verify a ZeroTest proof using pairing-based verification
    static bool verify(const ZeroTestProof& proof, const Fr& omega, 
                      size_t l, const KZG::SetupParams& params);
    
    // Enhanced verification with additional cryptographic checks
    static bool verify_with_full_checks(const ZeroTestProof& proof, const Fr& omega, 
                                       size_t l, const KZG::SetupParams& params);

private:
    // For debugging/testing
    static bool verify_division(const vector<Fr>& dividend, 
                               const vector<Fr>& divisor,
                               const vector<Fr>& quotient);
};