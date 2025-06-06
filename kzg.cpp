#include "kzg.h"
#include <random>
#include <iostream>
#include <cstring>
#include <cassert>
#include <stdexcept>

// Setup phase - generates structured reference string
KZG::SetupParams KZG::Setup(size_t max_degree) {
    SetupParams params;
    params.max_degree = max_degree;
    
    Fr tau;
    tau.setByCSPRNG();
    
    G1 g1;
    const char* g1_seed = "BN_SNARK1_G1_GENERATOR";
    hashAndMapToG1(g1, g1_seed, strlen(g1_seed));
    
    G2 g2;
    const char* g2_seed = "BN_SNARK1_G2_GENERATOR";
    hashAndMapToG2(g2, g2_seed, strlen(g2_seed));
    
    params.g1_powers.resize(max_degree + 1);
    Fr tau_power = Fr(1);
    
    for (size_t i = 0; i <= max_degree; i++) {
        G1::mul(params.g1_powers[i], g1, tau_power);
        if (i < max_degree) {
            Fr::mul(tau_power, tau_power, tau);
        }
    }
    
    params.g2_powers.resize(max_degree + 1);
    tau_power = Fr(1);
    for (size_t i = 0; i <= max_degree; i++) {
        G2::mul(params.g2_powers[i], g2, tau_power);
        if (i < max_degree) {
            Fr::mul(tau_power, tau_power, tau);
        }
    }
    
    tau.clear();
    tau_power.clear();
    
    return params;
}

// Commit to polynomial coefficients
KZG::Commitment KZG::Commit(const std::vector<Fr>& coefficients, const SetupParams& params) {
    if (coefficients.size() > params.max_degree + 1) {
        throw std::invalid_argument("Polynomial degree exceeds setup parameters");
    }
    
    Commitment commitment;
    commitment.commit.clear();
    
    std::vector<G1> bases;
    std::vector<Fr> scalars;
    
    for (size_t i = 0; i < coefficients.size(); i++) {
        if (!coefficients[i].isZero()) {
            bases.push_back(params.g1_powers[i]);
            scalars.push_back(coefficients[i]);
        }
    }
    
    if (!bases.empty()) {
        G1::mulVec(commitment.commit, bases.data(), scalars.data(), bases.size());
    }
    
    return commitment;
}

// Create witness for evaluation at point z
KZG::Proof KZG::CreateWitness(const std::vector<Fr>& coefficients, const Fr& z, const SetupParams& params) {
    Proof proof;
    
    proof.evaluation = evaluate_polynomial(coefficients, z);
    
    std::vector<Fr> f_minus_fz = coefficients;
    if (!f_minus_fz.empty()) {
        Fr::sub(f_minus_fz[0], f_minus_fz[0], proof.evaluation);
    }
    
    std::vector<Fr> quotient = divide_polynomial_by_linear(f_minus_fz, z);
    
    proof.witness.clear();
    for (size_t i = 0; i < quotient.size() && i < params.g1_powers.size(); i++) {
        if (!quotient[i].isZero()) {
            G1 term;
            G1::mul(term, params.g1_powers[i], quotient[i]);
            G1::add(proof.witness, proof.witness, term);
        }
    }
    
    return proof;
}

// Verify evaluation proof using pairing
bool KZG::VerifyEval(const Commitment& commitment, const Fr& z, const Proof& proof, const SetupParams& params) {
    if (params.g1_powers.empty() || params.g2_powers.size() < 2) {
        return false;
    }
    
    G1 g1_eval;
    G1::mul(g1_eval, params.g1_powers[0], proof.evaluation);
    G1 left_g1;
    G1::sub(left_g1, commitment.commit, g1_eval);
    
    G2 g2_z;
    G2::mul(g2_z, params.g2_powers[0], z);
    G2 right_g2;
    G2::sub(right_g2, params.g2_powers[1], g2_z);
    
    GT left_pairing, right_pairing;
    pairing(left_pairing, left_g1, params.g2_powers[0]);
    pairing(right_pairing, proof.witness, right_g2);
    
    return left_pairing == right_pairing;
}

// Evaluate polynomial at point x using Horner's method
Fr KZG::evaluate_polynomial(const std::vector<Fr>& coefficients, const Fr& x) {
    if (coefficients.empty()) {
        return Fr(0);
    }
    
    Fr result = coefficients.back();
    
    for (int i = coefficients.size() - 2; i >= 0; i--) {
        Fr::mul(result, result, x);
        Fr::add(result, result, coefficients[i]);
    }
    
    return result;
}

// Optimized division by linear polynomial (x - z)
std::vector<Fr> KZG::divide_polynomial_by_linear(const std::vector<Fr>& dividend, const Fr& z) {
    if (dividend.empty()) {
        return {};
    }
    
    if (dividend.size() == 1) {
        return {};
    }
    
    std::vector<Fr> quotient(dividend.size() - 1);
    
    quotient[quotient.size() - 1] = dividend.back();
    
    for (int i = quotient.size() - 2; i >= 0; i--) {
        Fr temp;
        Fr::mul(temp, quotient[i + 1], z);
        Fr::add(quotient[i], dividend[i + 1], temp);
    }
    
    return quotient;
}

// General polynomial division
std::vector<Fr> KZG::divide_polynomial(const std::vector<Fr>& dividend, const std::vector<Fr>& divisor) {
    if (divisor.empty()) {
        throw std::invalid_argument("Divisor cannot be empty");
    }
    
    std::vector<Fr> clean_divisor = divisor;
    while (!clean_divisor.empty() && clean_divisor.back().isZero()) {
        clean_divisor.pop_back();
    }
    
    if (clean_divisor.empty() || clean_divisor.back().isZero()) {
        throw std::invalid_argument("Divisor cannot be zero polynomial");
    }
    
    std::vector<Fr> clean_dividend = dividend;
    while (!clean_dividend.empty() && clean_dividend.back().isZero()) {
        clean_dividend.pop_back();
    }
    
    if (clean_dividend.empty()) {
        return {};
    }
    
    if (clean_dividend.size() < clean_divisor.size()) {
        return {};
    }
    
    size_t quotient_degree = clean_dividend.size() - clean_divisor.size();
    std::vector<Fr> quotient(quotient_degree + 1, Fr(0));
    std::vector<Fr> remainder = clean_dividend;
    
    Fr leading_coeff_inv;
    Fr::inv(leading_coeff_inv, clean_divisor.back());
    
    for (int i = quotient_degree; i >= 0; i--) {
        if (remainder.size() >= clean_divisor.size()) {
            Fr::mul(quotient[i], remainder.back(), leading_coeff_inv);
            
            for (size_t j = 0; j < clean_divisor.size(); j++) {
                if (remainder.size() >= clean_divisor.size() - j) {
                    size_t remainder_idx = remainder.size() - clean_divisor.size() + j;
                    Fr term;
                    Fr::mul(term, quotient[i], clean_divisor[j]);
                    Fr::sub(remainder[remainder_idx], remainder[remainder_idx], term);
                }
            }
            
            if (!remainder.empty() && remainder.back().isZero()) {
                remainder.pop_back();
            }
        }
    }
    
    return quotient;
}