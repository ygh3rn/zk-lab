#include "kzg.h"
#include <random>
#include <iostream>
#include <cassert>

// Helper function for G2 commitments (defined once here, used by both KZG and PIOP)
G2 commit_g2_internal(const Polynomial& poly, const std::vector<G2>& g2_powers) {
    const auto& coeffs = poly.get_coefficients();
    
    G2 commitment;
    commitment.clear(); // Initialize to zero
    
    // Compute commitment as sum of coeff_i * tau^i * g2
    for (size_t i = 0; i < coeffs.size() && i < g2_powers.size(); ++i) {
        if (!coeffs[i].isZero()) {
            G2 temp;
            G2::mul(temp, g2_powers[i], coeffs[i]);
            commitment += temp;
        }
    }
    
    return commitment;
}

KZG::KZG(size_t degree) : max_degree(degree) {
    // Initialize pairing for BN curve
    mcl::bn::initPairing();
    
    setup(degree);
}

void KZG::setup(size_t degree) {
    max_degree = degree;
    
    // Generate random tau (in practice, this should be done through a trusted setup)
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Generate random tau
    tau.setByCSPRNG();
    
    // Resize power arrays - we need more G2 powers for PIOP verification
    g1_powers.resize(degree + 1);
    
    // CRITICAL: Extend G2 powers beyond just [g2, tau*g2] for pairing verification
    // We need at least degree+1 powers for vanishing polynomial commitments
    size_t g2_size = std::max(size_t(2), degree + 1);
    g2_powers.resize(g2_size);
    
    // Set generators
    G1 g1;
    G2 g2;
    mapToG1(g1, 1); // Use deterministic generator
    mapToG2(g2, 1); // Use deterministic generator
    
    // Compute powers of tau in G1: [g1, tau*g1, tau^2*g1, ..., tau^degree*g1]
    g1_powers[0] = g1;
    Fr tau_power = Fr(1);
    
    for (size_t i = 1; i <= degree; ++i) {
        tau_power *= tau;
        G1::mul(g1_powers[i], g1, tau_power);
    }
    
    // Compute powers of tau in G2: [g2, tau*g2, tau^2*g2, ..., tau^degree*g2]
    // This is CRITICAL for pairing-based PIOP verification
    g2_powers[0] = g2;
    tau_power = Fr(1);
    
    for (size_t i = 1; i < g2_size; ++i) {
        tau_power *= tau;
        G2::mul(g2_powers[i], g2, tau_power);
    }
    
    std::cout << "KZG setup completed for degree " << degree << std::endl;
    std::cout << "  Generated " << g1_powers.size() << " G1 powers: [g₁, τg₁, τ²g₁, ..., τ^" << degree << "g₁]" << std::endl;
    std::cout << "  Generated " << g2_powers.size() << " G2 powers: [g₂, τg₂, τ²g₂, ..., τ^" << (g2_size-1) << "g₂]" << std::endl;
    std::cout << "  Ready for full pairing-based PIOP verification" << std::endl;
}

G1 KZG::commit(const Polynomial& poly) {
    const auto& coeffs = poly.get_coefficients();
    assert(coeffs.size() <= max_degree + 1 && "Polynomial degree exceeds setup degree");
    
    G1 commitment;
    commitment.clear(); // Initialize to zero
    
    // Compute commitment as sum of coeff_i * tau^i * g1
    for (size_t i = 0; i < coeffs.size(); ++i) {
        if (!coeffs[i].isZero()) {
            G1 temp;
            G1::mul(temp, g1_powers[i], coeffs[i]);
            commitment += temp;
        }
    }
    
    return commitment;
}

G1 KZG::create_witness(const Polynomial& poly, const Fr& point) {
    // Compute quotient polynomial q(x) = (f(x) - f(point)) / (x - point)
    Fr eval_at_point = poly.evaluate(point);
    
    // Create polynomial f(x) - f(point)
    Polynomial shifted_poly = poly - Polynomial({eval_at_point});
    
    // Use the optimized division for linear factors
    Polynomial quotient = shifted_poly.divide_by_linear(point);
    
    // Commit to quotient polynomial
    return commit(quotient);
}

bool KZG::verify_eval(const G1& commitment, const Fr& point, const Fr& value, const G1& witness) {
    // FULL PAIRING-BASED VERIFICATION as required by KZG paper
    // Verify that e(commitment - value*g1, g2) = e(witness, tau*g2 - point*g2)
    
    std::cout << "  KZG: Performing pairing-based evaluation verification" << std::endl;
    
    // Check if we have enough G2 powers
    if (g2_powers.size() < 2) {
        std::cout << "    ERROR: Not enough G2 powers for verification" << std::endl;
        return false;
    }
    
    // Compute commitment - value*g1
    G1 left_g1 = commitment;
    G1 value_g1;
    G1::mul(value_g1, g1_powers[0], value);
    left_g1 -= value_g1;
    
    // Compute tau*g2 - point*g2 = (tau - point)*g2
    G2 right_g2 = g2_powers[1]; // tau*g2
    G2 point_g2;
    G2::mul(point_g2, g2_powers[0], point);
    right_g2 -= point_g2;
    
    // Compute pairings
    Fp12 left_pairing, right_pairing;
    pairing(left_pairing, left_g1, g2_powers[0]);
    pairing(right_pairing, witness, right_g2);
    
    // Check if pairings are equal
    bool result = (left_pairing == right_pairing);
    
    if (result) {
        std::cout << "    ✓ KZG pairing verification PASSED" << std::endl;
    } else {
        std::cout << "    ✗ KZG pairing verification FAILED" << std::endl;
    }
    
    return result;
}

G1 KZG::create_batch_witness(const Polynomial& poly, const std::vector<Fr>& points) {
    // Compute interpolation polynomial for evaluations at given points
    std::vector<std::pair<Fr, Fr>> eval_points;
    eval_points.reserve(points.size());
    
    for (const Fr& point : points) {
        Fr value = poly.evaluate(point);
        eval_points.emplace_back(point, value);
    }
    
    // Create interpolation polynomial r(x)
    Polynomial r_poly = Polynomial::lagrange_interpolation(eval_points);
    
    // Create vanishing polynomial Z(x) = ∏(x - point_i)
    Polynomial vanishing_poly = Polynomial::from_roots(points);
    
    // Compute quotient (f(x) - r(x)) / Z(x) using the division operator
    Polynomial numerator = poly - r_poly;
    Polynomial quotient = numerator / vanishing_poly;
    
    return commit(quotient);
}

bool KZG::verify_batch_eval(const G1& commitment, const std::vector<Fr>& points, 
                           const std::vector<Fr>& values, const G1& witness) {
    assert(points.size() == values.size() && "Points and values must have same size");
    
    // Create interpolation polynomial r(x) from points and values
    std::vector<std::pair<Fr, Fr>> eval_points;
    eval_points.reserve(points.size());
    
    for (size_t i = 0; i < points.size(); ++i) {
        eval_points.emplace_back(points[i], values[i]);
    }
    
    Polynomial r_poly = Polynomial::lagrange_interpolation(eval_points);
    G1 r_commitment = commit(r_poly);
    
    // Create vanishing polynomial Z(x) = ∏(x - point_i)
    Polynomial vanishing_poly = Polynomial::from_roots(points);
    
    // Use the global G2 commitment function
    G2 vanishing_commitment_g2 = commit_g2_internal(vanishing_poly, g2_powers);
    
    // FULL PAIRING VERIFICATION: e(commitment - r_commitment, g2) = e(witness, vanishing_commitment_g2)
    std::cout << "  KZG: Performing batch pairing verification" << std::endl;
    
    G1 left_g1 = commitment;
    left_g1 -= r_commitment;
    
    Fp12 left_pairing, right_pairing;
    pairing(left_pairing, left_g1, g2_powers[0]);
    pairing(right_pairing, witness, vanishing_commitment_g2);
    
    bool result = (left_pairing == right_pairing);
    
    if (result) {
        std::cout << "    ✓ KZG batch verification PASSED" << std::endl;
    } else {
        std::cout << "    ✗ KZG batch verification FAILED" << std::endl;
    }
    
    return result;
}