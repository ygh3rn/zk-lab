#include "zerotest.h"
#include "polynomial.h"
#include <stdexcept>

ZeroTestProof ZeroTest::prove(const std::vector<Fr>& polynomial, const Fr& omega, 
                             size_t l, const KZG::SetupParams& params) {
    ZeroTestProof proof;
    proof.commitment = KZG::Commit(polynomial, params);
    
    Fr omega_power = Fr(1);
    for (size_t i = 0; i < l; i++) {
        Fr eval = Polynomial::evaluate(polynomial, omega_power);
        if (!eval.isZero()) {
            throw std::invalid_argument("Polynomial does not vanish on subgroup");
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    std::vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    std::vector<Fr> quotient = Polynomial::divide(polynomial, vanishing_poly);
    
    KZG::Commitment quotient_commit = KZG::Commit(quotient, params);
    proof.quotient_proof.witness = quotient_commit.commit;
    proof.quotient_proof.evaluation = Fr(0);
    
    return proof;
}

bool ZeroTest::verify(const ZeroTestProof& proof, const Fr& omega, 
                     size_t l, const KZG::SetupParams& params) {
    if (proof.commitment.commit.isZero()) {
        return true;
    }
    
    if (l == 0 || l >= params.g2_powers.size()) {
        return false;
    }
    
    GT left_pairing, right_pairing;
    pairing(left_pairing, proof.commitment.commit, params.g2_powers[0]);
    
    G2 vanishing_g2;
    G2::sub(vanishing_g2, params.g2_powers[l], params.g2_powers[0]);
    
    pairing(right_pairing, proof.quotient_proof.witness, vanishing_g2);
    
    return left_pairing == right_pairing;
}