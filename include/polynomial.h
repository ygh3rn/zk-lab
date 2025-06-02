#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <mcl/bn.hpp>
#include "ntt.h"

using namespace mcl::bn;

class Polynomial {
private:
    std::vector<Fr> coefficients;
    
public:
    Polynomial() = default;
    Polynomial(const std::vector<Fr>& coeffs);
    
    // Basic arithmetic operations
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;  // Using NTT
    Polynomial operator/(const Polynomial& divisor) const; // ADDED: Division operator
    
    // ADDED: Division operations with quotient and remainder
    std::pair<Polynomial, Polynomial> divmod(const Polynomial& divisor) const;
    
    // ADDED: Optimized special case divisions
    Polynomial divide_by_linear(const Fr& root) const;  // O(d) for (x - root)
    Polynomial divide_by_vanishing(size_t subgroup_size) const;  // O(d) for x^n - 1
    
    // Evaluation and interpolation
    Fr evaluate(const Fr& x) const;
    static Polynomial lagrange_interpolation(const std::vector<std::pair<Fr, Fr>>& points);
    
    // ADDED: Utility functions for cryptographic operations
    static Polynomial vanishing_polynomial(size_t subgroup_size);  // Z_H(x) = x^n - 1
    static Polynomial from_roots(const std::vector<Fr>& roots);    // ∏(x - root_i)
    
    // Getters/Setters
    const std::vector<Fr>& get_coefficients() const { return coefficients; }
    void set_coefficients(const std::vector<Fr>& coeffs) { coefficients = coeffs; }
    size_t degree() const;
    bool is_zero() const;
    
    // Utility functions
    void resize(size_t new_size);
    void print() const;
    
private:
    // ADDED: Internal division implementation methods
    std::pair<Polynomial, Polynomial> long_division(const Polynomial& divisor) const;
    void remove_leading_zeros();
    Polynomial divide_by_linear_synthetic(const Fr& root) const;
    Polynomial divide_general(const Polynomial& divisor) const;
};

#endif