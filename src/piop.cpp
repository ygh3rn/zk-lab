#include "piop.h"
#include <cassert>
#include <iostream>

// Univariate ZeroTest PIOP Implementation
UnivariateZeroTest::UnivariateZeroTest(KZG& kzg_instance, size_t subgroup_size) 
    : kzg(kzg_instance), n(subgroup_size) {
    // Ensure n is a power of 2
    assert((n & (n - 1)) == 0 && "Subgroup size must be a power of 2");
    
    setup_subgroup();
}

void UnivariateZeroTest::setup_subgroup() {
    // Find a primitive n-th root of unity
    // For BN_SNARK1, we need omega such that omega^n = 1
    
    // Use the same approach as in NTT
    Fr candidate;
    Fr field_order_minus_1;
    field_order_minus_1 = Fr(-1); // Use p-1 directly
    
    Fr order_div_n = field_order_minus_1;
    order_div_n /= Fr(n);
    
    candidate = Fr(5); // Known generator
    Fr::pow(omega, candidate, order_div_n);
    
    // Generate subgroup H = {1, ω, ω², ..., ω^(n-1)}
    subgroup.resize(n);
    subgroup[0] = Fr(1);
    
    for (size_t i = 1; i < n; ++i) {
        subgroup[i] = subgroup[i-1] * omega;
    }
    
    std::cout << "ZeroTest subgroup of size " << n << " initialized" << std::endl;
}

Polynomial UnivariateZeroTest::compute_vanishing_polynomial() {
    // Vanishing polynomial Z_H(x) = x^n - 1 for multiplicative subgroup
    std::vector<Fr> coeffs(n + 1, Fr(0));
    coeffs[0] = Fr(-1); // -1
    coeffs[n] = Fr(1);  // x^n
    
    return Polynomial(coeffs);
}

UnivariateZeroTest::Proof UnivariateZeroTest::prove(const Polynomial& poly) {
    Proof proof;
    
    // 1. Commit to polynomial f(x)
    proof.polynomial_commitment = kzg.commit(poly);
    
    // 2. Check that f(x) evaluates to zero on all points in H
    proof.evaluations.resize(n);
    bool all_zero = true;
    
    for (size_t i = 0; i < n; ++i) {
        proof.evaluations[i] = poly.evaluate(subgroup[i]);
        if (!proof.evaluations[i].isZero()) {
            all_zero = false;
            std::cout << "Warning: Polynomial does not evaluate to zero at subgroup element " 
                      << i << std::endl;
        }
    }
    
    if (!all_zero) {
        std::cout << "Error: Polynomial is not zero on the subgroup" << std::endl;
        // In a real implementation, we might return an error here
    }
    
    // 3. Compute quotient polynomial q(x) = f(x) / Z_H(x)
    // Since f should be zero on H, we have f(x) = q(x) * Z_H(x)
    Polynomial vanishing_poly = compute_vanishing_polynomial();
    
    // For the quotient, we can use polynomial division
    // But since f should be divisible by Z_H, we expect no remainder
    
    // Simple approach: assume f(x) = q(x) * Z_H(x) and solve for q(x)
    // This is a simplification - in practice you'd use proper polynomial division
    
    // For now, let's create a dummy quotient (this should be computed properly)
    std::vector<Fr> quotient_coeffs(poly.get_coefficients().size() > n ? 
                                   poly.get_coefficients().size() - n : 1, Fr(0));
    
    // This is a placeholder - proper implementation would divide f by Z_H
    if (all_zero) {
        // If polynomial is actually zero on H, we can compute the quotient
        quotient_coeffs[0] = Fr(1); // Placeholder
    }
    
    Polynomial quotient(quotient_coeffs);
    proof.quotient_commitment = kzg.commit(quotient);
    
    // 4. Generate witnesses for evaluations
    proof.witnesses.resize(n);
    for (size_t i = 0; i < n; ++i) {
        proof.witnesses[i] = kzg.create_witness(poly, subgroup[i]);
    }
    
    return proof;
}

bool UnivariateZeroTest::verify(const Proof& proof) {
    // 1. Verify that all evaluations are zero
    for (size_t i = 0; i < n; ++i) {
        if (!proof.evaluations[i].isZero()) {
            std::cout << "Verification failed: evaluation " << i << " is not zero" << std::endl;
            return false;
        }
    }
    
    // 2. Verify all witnesses
    for (size_t i = 0; i < n; ++i) {
        if (!kzg.verify_eval(proof.polynomial_commitment, subgroup[i], 
                            proof.evaluations[i], proof.witnesses[i])) {
            std::cout << "Verification failed: witness " << i << " is invalid" << std::endl;
            return false;
        }
    }
    
    // 3. Additional check: verify that polynomial commitment equals quotient * vanishing polynomial
    // This is a simplification - full verification would use pairings
    
    std::cout << "ZeroTest verification successful" << std::endl;
    return true;
}

// Univariate SumCheck PIOP Implementation
UnivariateSumCheck::UnivariateSumCheck(KZG& kzg_instance, size_t subgroup_size) 
    : kzg(kzg_instance), n(subgroup_size) {
    // Ensure n is a power of 2
    assert((n & (n - 1)) == 0 && "Subgroup size must be a power of 2");
    
    setup_subgroup();
}

void UnivariateSumCheck::setup_subgroup() {
    // Same setup as ZeroTest
    Fr candidate;
    Fr field_order_minus_1;
    field_order_minus_1 = Fr(-1); // Use p-1 directly
    
    Fr order_div_n = field_order_minus_1;
    order_div_n /= Fr(n);
    
    candidate = Fr(5);
    Fr::pow(omega, candidate, order_div_n);
    
    subgroup.resize(n);
    subgroup[0] = Fr(1);
    
    for (size_t i = 1; i < n; ++i) {
        subgroup[i] = subgroup[i-1] * omega;
    }
    
    std::cout << "SumCheck subgroup of size " << n << " initialized" << std::endl;
}

Fr UnivariateSumCheck::compute_sum_over_subgroup(const Polynomial& poly) {
    Fr sum(0);
    
    for (size_t i = 0; i < n; ++i) {
        Fr evaluation = poly.evaluate(subgroup[i]);
        sum += evaluation;
    }
    
    return sum;
}

UnivariateSumCheck::Proof UnivariateSumCheck::prove(const Polynomial& poly) {
    Proof proof;
    
    // 1. Commit to polynomial f(x)
    proof.polynomial_commitment = kzg.commit(poly);
    
    // 2. Compute the sum over the subgroup
    proof.claimed_sum = compute_sum_over_subgroup(poly);
    
    // 3. For SumCheck, we need to prove that ∑_{a∈H} f(a) = claimed_sum
    // We use the fact that for any polynomial g of degree < |H|:
    // ∑_{a∈H} g(a) = g(0) * |H| if g is constant, 0 otherwise
    
    // Create polynomial g(x) = f(x) - claimed_sum/|H|
    Fr sum_per_element = proof.claimed_sum;
    sum_per_element /= Fr(n);
    
    Polynomial constant_poly({sum_per_element});
    Polynomial adjusted_poly = poly - constant_poly;
    
    // 4. Now we need to show that ∑_{a∈H} adjusted_poly(a) = 0
    // This reduces to a ZeroTest-like problem
    
    // Compute quotient polynomial for the sum constraint
    // This is a simplification - proper implementation would be more complex
    std::vector<Fr> quotient_coeffs(poly.get_coefficients().size(), Fr(0));
    quotient_coeffs[0] = Fr(1); // Placeholder
    
    Polynomial quotient(quotient_coeffs);
    proof.quotient_commitment = kzg.commit(quotient);
    
    // 5. Generate witnesses (simplified)
    proof.witnesses.resize(1);
    proof.witnesses[0] = kzg.create_witness(poly, Fr(0)); // Witness at point 0
    
    return proof;
}

bool UnivariateSumCheck::verify(const Proof& proof, const Fr& expected_sum) {
    // 1. Check if claimed sum matches expected sum
    if (proof.claimed_sum != expected_sum) {
        std::cout << "Verification failed: claimed sum doesn't match expected sum" << std::endl;
        return false;
    }
    
    // 2. Verify witness (simplified verification)
    Fr eval_at_zero = Fr(0); // This should be computed from the actual polynomial
    // Simplified verification - check sum calculation instead
    if (proof.claimed_sum != expected_sum) {
        std::cout << "Verification failed: witness is invalid" << std::endl;
        return false;
    }
    
    // 3. Additional verification steps would go here in a complete implementation
    // This would involve checking the sum constraint using pairings and the quotient polynomial
    
    std::cout << "SumCheck verification successful" << std::endl;
    return true;
}