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
            throw std::invalid_argument("Polynomial does not vanish at ω^" + std::to_string(i));
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    std::vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    std::vector<Fr> quotient = Polynomial::divide(polynomial, vanishing_poly);
    
    std::vector<Fr> product = Polynomial::multiply(quotient, vanishing_poly);
    
    if (product.size() != polynomial.size()) {
        if (product.size() < polynomial.size()) {
            product.resize(polynomial.size(), Fr(0));
        } else {
            for (size_t i = polynomial.size(); i < product.size(); i++) {
                if (!product[i].isZero()) {
                    throw std::runtime_error("Polynomial division has non-zero remainder");
                }
            }
        }
    }
    
    for (size_t i = 0; i < polynomial.size(); i++) {
        if (!(polynomial[i] == product[i])) {
            throw std::runtime_error("Polynomial division verification failed");
        }
    }
    
    KZG::Commitment quotient_commit = KZG::Commit(quotient, params);
    proof.quotient_proof.witness = quotient_commit.commit;
    proof.quotient_proof.evaluation = Fr(0);
    
    return proof;
}

bool ZeroTest::verify(const ZeroTestProof& proof, const Fr& omega, 
                     size_t l, const KZG::SetupParams& params) {
    if (l == 0 || l >= params.g2_powers.size()) {
        return false;
    }
    
    if (proof.commitment.commit.isZero()) {
        return proof.quotient_proof.witness.isZero();
    }
    
    if (l >= params.g2_powers.size()) {
        return false;
    }
    
    G2 vanishing_g2;
    G2::sub(vanishing_g2, params.g2_powers[l], params.g2_powers[0]);
    
    GT left_pairing, right_pairing;
    pairing(left_pairing, proof.commitment.commit, params.g2_powers[0]);
    pairing(right_pairing, proof.quotient_proof.witness, vanishing_g2);
    
    return left_pairing == right_pairing;
}

bool ZeroTest::verify_with_full_checks(const ZeroTestProof& proof, const Fr& omega, 
                                      size_t l, const KZG::SetupParams& params) {
    if (!verify(proof, omega, l, params)) {
        return false;
    }
    
    Fr omega_to_l;
    Fr::pow(omega_to_l, omega, Fr(l));
    if (!(omega_to_l == Fr(1))) {
        return false;
    }
    
    Fr omega_power = omega;
    for (size_t k = 1; k < l; k++) {
        if (omega_power == Fr(1)) {
            return false;
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    if (l > params.max_degree) {
        return false;
    }
    
    if (!proof.commitment.commit.isValidOrder()) {
        return false;
    }
    
    if (!proof.quotient_proof.witness.isValidOrder()) {
        return false;
    }
    
    return true;
}

// Helper function to verify polynomial division (for testing/debugging)
bool ZeroTest::verify_division(const std::vector<Fr>& dividend, 
                              const std::vector<Fr>& divisor,
                              const std::vector<Fr>& quotient) {
    std::vector<Fr> product = Polynomial::multiply(quotient, divisor);
    
    if (product.size() != dividend.size()) {
        if (product.size() < dividend.size()) {
            product.resize(dividend.size(), Fr(0));
        } else {
            for (size_t i = dividend.size(); i < product.size(); i++) {
                if (!product[i].isZero()) {
                    return false;
                }
            }
        }
    }
    
    for (size_t i = 0; i < dividend.size(); i++) {
        if (!(dividend[i] == product[i])) {
            return false;
        }
    }
    
    return true;
}