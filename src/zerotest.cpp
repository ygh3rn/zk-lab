#include "zerotest.h"
#include "polynomial.h"
#include <stdexcept>

pair<KZG::Commitment, KZG::Commitment> ZeroTest::commit_phase(
    const vector<Fr>& polynomial, const Fr& omega, size_t l, 
    const KZG::SetupParams& params) {
    
    Fr omega_power = Fr(1);
    for (size_t i = 0; i < l; i++) {
        Fr eval = Polynomial::evaluate(polynomial, omega_power);
        if (!eval.isZero()) {
            throw invalid_argument("Polynomial does not vanish at ω^" + to_string(i));
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    vector<Fr> quotient = Polynomial::divide(polynomial, vanishing_poly);
    
    vector<Fr> product = Polynomial::multiply(quotient, vanishing_poly);
    if (product.size() != polynomial.size()) {
        throw runtime_error("Polynomial division has non-zero remainder");
    }
    for (size_t i = 0; i < polynomial.size(); i++) {
        if (!(polynomial[i] == product[i])) {
            throw runtime_error("Polynomial division has non-zero remainder");
        }
    }
    
    KZG::Commitment f_commit = KZG::Commit(polynomial, params);
    KZG::Commitment q_commit = KZG::Commit(quotient, params);
    
    return {f_commit, q_commit};
}

ZeroTestProof ZeroTest::prove(const vector<Fr>& polynomial, 
                             const Fr& challenge, const Fr& omega, size_t l,
                             const KZG::SetupParams& params) {
    ZeroTestProof proof;
    
    auto [f_commit, q_commit] = commit_phase(polynomial, omega, l, params);
    proof.f_commitment = f_commit;
    proof.q_commitment = q_commit;
    proof.challenge = challenge;
    
    vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    vector<Fr> quotient = Polynomial::divide(polynomial, vanishing_poly);
    
    proof.f_eval = Polynomial::evaluate(polynomial, challenge);
    proof.q_eval = Polynomial::evaluate(quotient, challenge);
    
    proof.f_eval_proof = KZG::CreateWitness(polynomial, challenge, params);
    proof.q_eval_proof = KZG::CreateWitness(quotient, challenge, params);
    
    return proof;
}

bool ZeroTest::verify(const ZeroTestProof& proof, const Fr& omega, size_t l,
                     const KZG::SetupParams& params) {
    if (!KZG::VerifyEval(proof.f_commitment, proof.challenge, proof.f_eval_proof, params)) {
        return false;
    }
    
    if (!KZG::VerifyEval(proof.q_commitment, proof.challenge, proof.q_eval_proof, params)) {
        return false;
    }
    
    Fr vanishing_eval;
    Fr::pow(vanishing_eval, proof.challenge, Fr(l));
    Fr::sub(vanishing_eval, vanishing_eval, Fr(1));
    
    Fr expected_f_eval;
    Fr::mul(expected_f_eval, proof.q_eval, vanishing_eval);
    
    return proof.f_eval == expected_f_eval;
}

bool ZeroTest::verify_with_full_checks(const ZeroTestProof& proof, const Fr& omega, 
                                      size_t l, const KZG::SetupParams& params) {
    if (!verify(proof, omega, l, params)) {
        return false;
    }
    
    Fr omega_to_l;
    Fr::pow(omega_to_l, omega, Fr(l));
    if (omega_to_l != Fr(1)) {
        return false;
    }
    
    Fr omega_power = omega;
    for (size_t k = 1; k < l; k++) {
        if (omega_power == Fr(1)) {
            return false;
        }
        Fr::mul(omega_power, omega_power, omega);
    }
    
    if (!proof.f_commitment.commit.isValidOrder() || 
        !proof.q_commitment.commit.isValidOrder()) {
        return false;
    }
    
    return true;
}