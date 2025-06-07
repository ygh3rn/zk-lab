#include "sumcheck.h"
#include "polynomial.h"

SumCheckProof SumCheck::prove(const std::vector<Fr>& polynomial, const Fr& omega, 
                             size_t l, const KZG::SetupParams& params) {
    SumCheckProof proof;
    proof.commitment = KZG::Commit(polynomial, params);
    
    // Compute actual sum over subgroup
    proof.claimed_sum = Polynomial::sum_on_subgroup(polynomial, omega, l);
    
    // Create adjusted polynomial: f_adj(x) = f(x) - (sum/l)
    // This ensures sum of f_adj over subgroup is 0
    std::vector<Fr> adjusted_poly = polynomial;
    
    Fr avg_value;
    Fr::div(avg_value, proof.claimed_sum, Fr(l));
    
    if (adjusted_poly.empty()) {
        adjusted_poly = {Fr(0)};
    }
    Fr::sub(adjusted_poly[0], adjusted_poly[0], avg_value);
    
    // Compute quotient: q(x) = f_adj(x) / Z_H(x)
    std::vector<Fr> vanishing_poly = Polynomial::vanishing(l);
    std::vector<Fr> quotient = Polynomial::divide(adjusted_poly, vanishing_poly);
    
    // Commit to quotient
    KZG::Commitment quotient_commit = KZG::Commit(quotient, params);
    proof.adjusted_proof.witness = quotient_commit.commit;
    proof.adjusted_proof.evaluation = Fr(0);
    
    return proof;
}

bool SumCheck::verify(const SumCheckProof& proof, const Fr& omega, 
                     size_t l, const KZG::SetupParams& params) {
    if (l == 0 || l >= params.g2_powers.size()) {
        return false;
    }
    
    // Compute adjusted commitment: C_adj = C - (sum/l) * g₁
    Fr avg_value;
    Fr::div(avg_value, proof.claimed_sum, Fr(l));
    
    G1 const_commit;
    G1::mul(const_commit, params.g1_powers[0], avg_value);
    
    G1 adjusted_commit;
    G1::sub(adjusted_commit, proof.commitment.commit, const_commit);
    
    // Pairing check: e(C_adj, g₂) = e(q(τ), Z_H(τ)g₂)
    GT left_pairing, right_pairing;
    pairing(left_pairing, adjusted_commit, params.g2_powers[0]);
    
    // Compute Z_H(τ) = τ^l - 1 in G2
    G2 vanishing_g2;
    G2::sub(vanishing_g2, params.g2_powers[l], params.g2_powers[0]);
    
    pairing(right_pairing, proof.adjusted_proof.witness, vanishing_g2);
    
    return left_pairing == right_pairing;
}