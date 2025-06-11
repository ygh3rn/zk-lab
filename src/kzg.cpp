#include "kzg.h"
#include "polynomial.h"
#include <stdexcept>
#include <cstring>

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
    params.g2_powers.resize(max_degree + 1);
    
    Fr tau_power = Fr(1);
    for (size_t i = 0; i <= max_degree; i++) {
        G1::mul(params.g1_powers[i], g1, tau_power);
        G2::mul(params.g2_powers[i], g2, tau_power);
        if (i < max_degree) {
            Fr::mul(tau_power, tau_power, tau);
        }
    }
    
    tau.clear();
    
    return params;
}

KZG::Commitment KZG::Commit(const vector<Fr>& coefficients, const SetupParams& params) {
    if (coefficients.size() > params.max_degree + 1) {
        throw invalid_argument("Polynomial degree exceeds setup parameters");
    }
    
    Commitment commitment;
    commitment.commit.clear();
    
    vector<G1> bases;
    vector<Fr> scalars;
    
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

KZG::Proof KZG::CreateWitness(const vector<Fr>& coefficients, const Fr& z, const SetupParams& params) {
    Proof proof;
    
    proof.evaluation = Polynomial::evaluate(coefficients, z);
    
    vector<Fr> f_minus_fz = coefficients;
    if (!f_minus_fz.empty()) {
        Fr::sub(f_minus_fz[0], f_minus_fz[0], proof.evaluation);
    }
    
    vector<Fr> quotient = Polynomial::divide_by_linear(f_minus_fz, z);
    
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

KZG::BatchProof KZG::CreateWitnessBatch(const vector<Fr>& coefficients, 
                                        const vector<Fr>& eval_points, 
                                        const SetupParams& params) {
    if (eval_points.empty()) {
        throw invalid_argument("Evaluation points cannot be empty");
    }
    
    if (coefficients.size() > params.max_degree + 1) {
        throw invalid_argument("Polynomial degree exceeds setup parameters");
    }
    
    BatchProof proof;
    proof.evaluation_points = eval_points;
    proof.evaluations.resize(eval_points.size());
    
    for (size_t i = 0; i < eval_points.size(); i++) {
        proof.evaluations[i] = Polynomial::evaluate(coefficients, eval_points[i]);
    }
    
    vector<Fr> remainder_poly = Polynomial::interpolate_lagrange(eval_points, proof.evaluations);
    
    vector<Fr> f_minus_r = coefficients;
    size_t max_size = max(f_minus_r.size(), remainder_poly.size());
    f_minus_r.resize(max_size, Fr(0));
    remainder_poly.resize(max_size, Fr(0));
    
    for (size_t i = 0; i < max_size; i++) {
        Fr::sub(f_minus_r[i], f_minus_r[i], remainder_poly[i]);
    }
    
    vector<Fr> product_poly = {Fr(1)};
    for (const Fr& point : eval_points) {
        vector<Fr> linear = {Fr(0), Fr(1)};
        Fr::sub(linear[0], linear[0], point);
        product_poly = Polynomial::multiply(product_poly, linear);
    }
    
    vector<Fr> quotient = Polynomial::divide(f_minus_r, product_poly);
    
    vector<Fr> verification = Polynomial::multiply(quotient, product_poly);
    if (verification.size() != f_minus_r.size()) {
        verification.resize(f_minus_r.size(), Fr(0));
    }
    for (size_t i = 0; i < f_minus_r.size(); i++) {
        if (!(f_minus_r[i] == verification[i])) {
            throw runtime_error("Batch witness creation failed: division not exact");
        }
    }
    
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

bool KZG::VerifyEvalBatch(const Commitment& commitment, const BatchProof& proof, 
                         const SetupParams& params) {
    if (proof.evaluation_points.empty() || 
        proof.evaluation_points.size() != proof.evaluations.size()) {
        return false;
    }
    
    if (proof.evaluation_points.size() > params.max_degree) {
        return false;
    }
    
    vector<Fr> remainder_poly = Polynomial::interpolate_lagrange(
        proof.evaluation_points, proof.evaluations);
    
    G1 r_commit;
    r_commit.clear();
    for (size_t i = 0; i < remainder_poly.size() && i < params.g1_powers.size(); i++) {
        if (!remainder_poly[i].isZero()) {
            G1 term;
            G1::mul(term, params.g1_powers[i], remainder_poly[i]);
            G1::add(r_commit, r_commit, term);
        }
    }
    
    G1 left_g1;
    G1::sub(left_g1, commitment.commit, r_commit);
    
    vector<Fr> product_poly = {Fr(1)};
    for (const Fr& point : proof.evaluation_points) {
        vector<Fr> linear = {Fr(0), Fr(1)};
        Fr::sub(linear[0], linear[0], point);
        product_poly = Polynomial::multiply(product_poly, linear);
    }
    
    G2 product_g2;
    product_g2.clear();
    for (size_t i = 0; i < product_poly.size() && i < params.g2_powers.size(); i++) {
        if (!product_poly[i].isZero()) {
            G2 term;
            G2::mul(term, params.g2_powers[i], product_poly[i]);
            G2::add(product_g2, product_g2, term);
        }
    }
    
    GT left_pairing, right_pairing;
    pairing(left_pairing, left_g1, params.g2_powers[0]);
    pairing(right_pairing, proof.witness, product_g2);
    
    return left_pairing == right_pairing;
}