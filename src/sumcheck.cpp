#include "sumcheck.h"
#include "polynomial.h"
#include <stdexcept>

pair<KZG::Commitment, KZG::Commitment> SumCheck::commit_phase(
    const vector<Fr>& polynomial, const Fr& omega, size_t l, 
    const KZG::SetupParams& params) {
    
    Fr sum = Polynomial::sum_on_subgroup(polynomial, omega, l);
    if (!sum.isZero()) {
        throw invalid_argument("Polynomial sum over subgroup is not zero: " + sum.getStr());
    }
    
    if (!polynomial.empty() && !polynomial[0].isZero()) {
        throw invalid_argument("Polynomial has non-zero constant term, cannot have sum=0");
    }
    
    vector<Fr> quotient;
    if (polynomial.size() <= 1) {
        quotient = {Fr(0)};
    } else {
        quotient.resize(polynomial.size() - 1);
        for (size_t i = 1; i < polynomial.size(); i++) {
            quotient[i - 1] = polynomial[i];
        }
    }
    
    vector<Fr> product = Polynomial::multiply(quotient, {Fr(0), Fr(1)});
    if (product.size() < polynomial.size()) {
        product.resize(polynomial.size(), Fr(0));
    }
    
    bool division_valid = true;
    size_t check_size = max(product.size(), polynomial.size());
    for (size_t i = 0; i < check_size; i++) {
        Fr poly_coeff = (i < polynomial.size()) ? polynomial[i] : Fr(0);
        Fr prod_coeff = (i < product.size()) ? product[i] : Fr(0);
        if (!(poly_coeff == prod_coeff)) {
            division_valid = false;
            break;
        }
    }
    
    if (!division_valid) {
        throw runtime_error("Polynomial is not divisible by X");
    }
    
    KZG::Commitment f_commit = KZG::Commit(polynomial, params);
    KZG::Commitment q_commit = KZG::Commit(quotient, params);
    
    return {f_commit, q_commit};
}

SumCheckProof SumCheck::prove(const vector<Fr>& polynomial, 
                             const Fr& challenge, const Fr& omega, size_t l,
                             const KZG::SetupParams& params) {
    SumCheckProof proof;
    
    proof.claimed_sum = Polynomial::sum_on_subgroup(polynomial, omega, l);
    if (!proof.claimed_sum.isZero()) {
        throw invalid_argument("This implementation only supports sum = 0 case");
    }
    
    auto [f_commit, q_commit] = commit_phase(polynomial, omega, l, params);
    proof.f_commitment = f_commit;
    proof.q_commitment = q_commit;
    proof.challenge = challenge;
    
    vector<Fr> quotient;
    if (polynomial.size() <= 1) {
        quotient = {Fr(0)};
    } else {
        quotient.resize(polynomial.size() - 1);
        for (size_t i = 1; i < polynomial.size(); i++) {
            quotient[i - 1] = polynomial[i];
        }
    }
    
    proof.f_eval = Polynomial::evaluate(polynomial, challenge);
    proof.q_eval = Polynomial::evaluate(quotient, challenge);
    
    proof.f_eval_proof = KZG::CreateWitness(polynomial, challenge, params);
    proof.q_eval_proof = KZG::CreateWitness(quotient, challenge, params);
    
    return proof;
}

bool SumCheck::verify(const SumCheckProof& proof, const Fr& omega, size_t l,
                     const KZG::SetupParams& params) {
    if (!proof.claimed_sum.isZero()) {
        return false;
    }
    
    if (!KZG::VerifyEval(proof.f_commitment, proof.challenge, proof.f_eval_proof, params)) {
        return false;
    }
    
    if (!KZG::VerifyEval(proof.q_commitment, proof.challenge, proof.q_eval_proof, params)) {
        return false;
    }
    
    Fr expected_f_eval;
    Fr::mul(expected_f_eval, proof.challenge, proof.q_eval);
    
    return proof.f_eval == expected_f_eval;
}

bool SumCheck::verify_with_full_checks(const SumCheckProof& proof, const Fr& omega, 
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
    
    if (l > params.max_degree) {
        return false;
    }
    
    if (!proof.f_commitment.commit.isValidOrder() || 
        !proof.q_commitment.commit.isValidOrder()) {
        return false;
    }
    
    return true;
}