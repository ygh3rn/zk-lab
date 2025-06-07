#include "zerotest.h"
#include "polynomial.h"
#include <stdexcept>

ZeroTestProof ZeroTest::prove(const std::vector<Fr>& polynomial, const Fr& omega, 
                             size_t l, const KZG::SetupParams& params) {
    ZeroTestProof proof;
    proof.commitment = KZG::Commit(polynomial, params);
    
    // Verify polynomial actually vanishes on subgroup
    Fr omega_power = Fr(1);
    for (size_t i = 0; i < l; i++) {
        Fr eval = Polynomial::evaluate(polynomial, omega_power);
        if (!eval.isZero()) {
            throw std::invalid_argument("Polynomial does not vanish on subgroup");
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    // Compute quotient polynomial q(x) = f(x) / Z_H(x)
    // where Z_H(x) = x^l - 1 is the vanishing polynomial
    std::vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    std::vector<Fr> quotient = Polynomial::divide(polynomial, vanishing_poly);
    
    // Commit to quotient
    KZG::Commitment quotient_commit = KZG::Commit(quotient, params);
    proof.quotient_proof.witness = quotient_commit.commit;
    proof.quotient_proof.evaluation = Fr(0);
    
    return proof;
}

bool ZeroTest::verify(const ZeroTestProof& proof, const Fr& omega, 
                     size_t l, const KZG::SetupParams& params) {
    if (proof.commitment.commit.isZero()) {
        return true; // Zero polynomial trivially vanishes
    }
    
    if (l == 0 || l >= params.g2_powers.size()) {
        return false;
    }
    
    // Pairing check: e(f(τ), g₂) = e(q(τ), Z_H(τ)g₂)
    // where Z_H(τ) = τ^l - 1
    GT left_pairing, right_pairing;
    pairing(left_pairing, proof.commitment.commit, params.g2_powers[0]);
    
    // Compute Z_H(τ) = τ^l - 1 in G2
    G2 vanishing_g2;
    G2::sub(vanishing_g2, params.g2_powers[l], params.g2_powers[0]);
    
    pairing(right_pairing, proof.quotient_proof.witness, vanishing_g2);
    
    return left_pairing == right_pairing;
}