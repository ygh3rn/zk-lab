#include "kzg.h"
#include <random>
#include <iostream>
#include <cstring>
#include <cassert>

// Setup phase - generates structured reference string
KZG::SetupParams KZG::Setup(size_t max_degree) {
    SetupParams params;
    params.max_degree = max_degree;
    
    // Generate random secret τ (toxic waste)
    Fr tau;
    tau.setByCSPRNG();
    
    // Generate G1 generator
    G1 g1;
    const char* g1_seed = "BN_SNARK1_G1_GENERATOR";
    hashAndMapToG1(g1, g1_seed, strlen(g1_seed));
    
    // Generate G2 generator  
    G2 g2;
    const char* g2_seed = "BN_SNARK1_G2_GENERATOR";
    hashAndMapToG2(g2, g2_seed, strlen(g2_seed));
    
    // Compute G1 powers: [g₁, g₁^τ, g₁^τ², ..., g₁^τⁿ]
    params.g1_powers.resize(max_degree + 1);
    Fr tau_power = Fr(1);  // τ⁰ = 1
    
    for (size_t i = 0; i <= max_degree; i++) {
        G1::mul(params.g1_powers[i], g1, tau_power);
        if (i < max_degree) {
            Fr::mul(tau_power, tau_power, tau);  // τⁱ⁺¹ = τⁱ * τ
        }
    }
    
    // Compute G2 powers: [g₂, g₂^τ] (only need 2 for verification)
    params.g2_powers.resize(2);
    G2::mul(params.g2_powers[0], g2, Fr(1));     // g₂
    G2::mul(params.g2_powers[1], g2, tau);       // g₂^τ
    
    // Toxic waste disposal - τ is automatically destroyed when leaving scope
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
    commitment.commit.clear();  // Initialize to identity
    
    // Compute C = Σᵢ(aᵢ * g₁^(τⁱ)) using multi-scalar multiplication
    std::vector<G1> bases;
    std::vector<Fr> scalars;
    
    for (size_t i = 0; i < coefficients.size(); i++) {
        if (!coefficients[i].isZero()) {  // Skip zero coefficients
            bases.push_back(params.g1_powers[i]);
            scalars.push_back(coefficients[i]);
        }
    }
    
    if (!bases.empty()) {
        // Use MCL's optimized multi-scalar multiplication
        G1::mulVec(commitment.commit, bases.data(), scalars.data(), bases.size());
    }
    
    return commitment;
}

// Create witness for evaluation at point z
KZG::Proof KZG::CreateWitness(const std::vector<Fr>& coefficients, const Fr& z, const SetupParams& params) {
    Proof proof;
    
    // Evaluate f(z) = Σᵢ(aᵢ * zⁱ)
    proof.evaluation = evaluate_polynomial(coefficients, z);
    
    // Compute quotient polynomial q(x) = (f(x) - f(z)) / (x - z)
    std::vector<Fr> f_minus_fz = coefficients;
    if (!f_minus_fz.empty()) {
        Fr::sub(f_minus_fz[0], f_minus_fz[0], proof.evaluation);
    }
    
    // Polynomial division by (x - z)
    std::vector<Fr> quotient = divide_polynomial_by_linear(f_minus_fz, z);
    
    // Commit to quotient polynomial: W = Σᵢ(qᵢ * g₁^(τⁱ))
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
    
    // Pairing verification: e(C - g₁^f(z), g₂) = e(W, g₂^τ - g₂^z)
    
    // Compute C - g₁^f(z)
    G1 g1_eval;
    G1::mul(g1_eval, params.g1_powers[0], proof.evaluation);
    G1 left_g1;
    G1::sub(left_g1, commitment.commit, g1_eval);
    
    // Compute g₂^τ - g₂^z  
    G2 g2_z;
    G2::mul(g2_z, params.g2_powers[0], z);
    G2 right_g2;
    G2::sub(right_g2, params.g2_powers[1], g2_z);
    
    // Compute pairings
    GT left_pairing, right_pairing;
    pairing(left_pairing, left_g1, params.g2_powers[0]);     // e(C - g₁^f(z), g₂)
    pairing(right_pairing, proof.witness, right_g2);         // e(W, g₂^τ - g₂^z)
    
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
        // Constant polynomial: should be zero if divisible
        return {};
    }
    
    std::vector<Fr> quotient(dividend.size() - 1);
    
    // Synthetic division algorithm
    quotient[quotient.size() - 1] = dividend.back();
    
    for (int i = quotient.size() - 2; i >= 0; i--) {
        Fr temp;
        Fr::mul(temp, quotient[i + 1], z);
        Fr::add(quotient[i], dividend[i + 1], temp);
    }
    
    return quotient;
}

// General polynomial division (kept for compatibility)
std::vector<Fr> KZG::divide_polynomial(const std::vector<Fr>& dividend, const std::vector<Fr>& divisor) {
    if (divisor.empty() || divisor.back().isZero()) {
        throw std::invalid_argument("Invalid divisor");
    }
    
    if (dividend.size() < divisor.size()) {
        return {};
    }
    
    std::vector<Fr> quotient(dividend.size() - divisor.size() + 1);
    std::vector<Fr> remainder = dividend;
    
    for (int i = quotient.size() - 1; i >= 0; i--) {
        if (remainder.size() >= divisor.size()) {
            Fr::div(quotient[i], remainder.back(), divisor.back());
            
            for (size_t j = 0; j < divisor.size(); j++) {
                Fr term;
                Fr::mul(term, quotient[i], divisor[j]);
                size_t rem_idx = remainder.size() - divisor.size() + j;
                Fr::sub(remainder[rem_idx], remainder[rem_idx], term);
            }
            
            remainder.pop_back();
        }
    }
    
    return quotient;
}