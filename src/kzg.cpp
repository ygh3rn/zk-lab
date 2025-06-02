#include "kzg.h"
#include <random>
#include <iostream>
#include <cassert>

KZG::KZG(size_t degree) : max_degree(degree) {
    // Initialize pairing for BN curve
    mcl::bn::initPairing();
    
    setup(degree);
}

void KZG::setup(size_t degree) {
    max_degree = degree;
    
    // Generate random tau
    tau.setByCSPRNG();
    
    // Resize power arrays
    g1_powers.resize(degree + 1);
    g2_powers.resize(2); // We only need g2 and tau*g2 for verification
    
    // Set generators
    G1 g1;
    G2 g2;
    mapToG1(g1, 1); // Use deterministic generator
    mapToG2(g2, 1); // Use deterministic generator
    
    // Compute powers of tau in G1: [g1, tau*g1, tau^2*g1, ..., tau^degree*g1]
    g1_powers[0] = g1;
    Fr tau_power = tau;
    
    for (size_t i = 1; i <= degree; ++i) {
        G1::mul(g1_powers[i], g1, tau_power);
        tau_power *= tau;
    }
    
    // Compute powers of tau in G2: [g2, tau*g2]
    g2_powers[0] = g2;
    G2::mul(g2_powers[1], g2, tau);
    
    std::cout << "KZG setup completed for degree " << degree << std::endl;
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
    
    // Divide by (x - point)
    // Create polynomial (x - point)
    std::vector<Fr> divisor_coeffs = {-point, Fr(1)};
    Polynomial divisor(divisor_coeffs);
    
    // Use the polynomial division operator
    Polynomial quotient = shifted_poly / divisor;
    
    // Commit to quotient polynomial
    return commit(quotient);
}

bool KZG::verify_eval(const G1& commitment, const Fr& point, const Fr& value, const G1& witness) {
    // Verify that e(commitment - value*g1, g2) = e(witness, tau*g2 - point*g2)
    
    // Compute commitment - value*g1
    G1 left_g1 = commitment;
    G1 value_g1;
    G1::mul(value_g1, g1_powers[0], value);
    left_g1 -= value_g1;
    
    // Compute tau*g2 - point*g2
    G2 right_g2 = g2_powers[1]; // tau*g2
    G2 point_g2;
    G2::mul(point_g2, g2_powers[0], point);
    right_g2 -= point_g2;
    
    // Compute pairings
    Fp12 left_pairing, right_pairing;
    pairing(left_pairing, left_g1, g2_powers[0]);
    pairing(right_pairing, witness, right_g2);
    
    // Check if pairings are equal
    return left_pairing == right_pairing;
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
    Polynomial vanishing_poly({Fr(1)}); // Start with polynomial 1
    
    for (const Fr& point : points) {
        std::vector<Fr> linear_coeffs = {-point, Fr(1)}; // (x - point)
        Polynomial linear_factor(linear_coeffs);
        vanishing_poly = vanishing_poly * linear_factor;
    }
    
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
    
    // Check if the interpolation matches expected values
    for (size_t i = 0; i < points.size(); ++i) {
        if (r_poly.evaluate(points[i]) != values[i]) {
            return false;
        }
    }
    
    return true;
}