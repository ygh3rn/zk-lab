#include "sumcheck.h"
#include "polynomial.h"
#include <stdexcept>

SumCheckProof SumCheck::prove(const vector<Fr>& polynomial, const Fr& omega, 
                             size_t l, const KZG::SetupParams& params) {
    SumCheckProof proof;
    proof.commitment = KZG::Commit(polynomial, params);
    
    proof.claimed_sum = Polynomial::sum_on_subgroup(polynomial, omega, l);
    
    vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    vector<Fr> quotient_div = Polynomial::divide(polynomial, vanishing_poly);
    
    vector<Fr> q_times_vanishing = Polynomial::multiply(quotient_div, vanishing_poly);
    
    vector<Fr> remainder = polynomial;
    if (remainder.size() < q_times_vanishing.size()) {
        remainder.resize(q_times_vanishing.size(), Fr(0));
    }
    
    for (size_t i = 0; i < q_times_vanishing.size(); i++) {
        Fr::sub(remainder[i], remainder[i], q_times_vanishing[i]);
    }
    
    while (!remainder.empty() && remainder.back().isZero()) {
        remainder.pop_back();
    }
    
    Fr r_constant = remainder.empty() ? Fr(0) : remainder[0];
    Fr expected_constant;
    Fr::div(expected_constant, proof.claimed_sum, Fr(l));
    
    if (!(r_constant == expected_constant)) {
        throw runtime_error("Invalid sum: computed sum doesn't match polynomial structure");
    }
    
    if (!proof.claimed_sum.isZero()) {
        throw invalid_argument("This implementation only supports sum = 0 case");
    }
    
    vector<Fr> quotient_x;
    if (remainder.empty() || remainder[0].isZero()) {
        if (remainder.size() <= 1) {
            quotient_x = {Fr(0)};
        } else {
            quotient_x.resize(remainder.size() - 1);
            for (size_t i = 1; i < remainder.size(); i++) {
                quotient_x[i - 1] = remainder[i];
            }
        }
    } else {
        throw runtime_error("Remainder has non-zero constant term, cannot divide by x");
    }
    
    KZG::Commitment quotient_commit = KZG::Commit(quotient_x, params);
    proof.adjusted_proof.witness = quotient_commit.commit;
    proof.adjusted_proof.evaluation = Fr(0);
    
    return proof;
}

bool SumCheck::verify(const SumCheckProof& proof, const Fr& omega, 
                     size_t l, const KZG::SetupParams& params) {
    if (l == 0 || params.g2_powers.size() < 2) {
        return false;
    }
    
    if (!proof.claimed_sum.isZero()) {
        return false;
    }
    
    if (proof.commitment.commit.isZero()) {
        return proof.adjusted_proof.witness.isZero();
    }
    
    GT left_pairing, right_pairing;
    pairing(left_pairing, proof.commitment.commit, params.g2_powers[0]);
    pairing(right_pairing, proof.adjusted_proof.witness, params.g2_powers[1]);
    
    bool pairing_check = (left_pairing == right_pairing);
    bool witness_valid = !proof.adjusted_proof.witness.isZero() || proof.commitment.commit.isZero();
    
    return pairing_check && witness_valid;
}

bool SumCheck::verify_with_full_checks(const SumCheckProof& proof, const Fr& omega, 
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
    
    if (!proof.adjusted_proof.witness.isValidOrder()) {
        return false;
    }
    
    return true;
}