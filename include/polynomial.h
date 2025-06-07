#pragma once

#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;

class Polynomial {
public:
    // Evaluate polynomial at point x using Horner's method
    static Fr evaluate(const std::vector<Fr>& coefficients, const Fr& x);
    
    // Multiply polynomials using NTT
    static std::vector<Fr> multiply(const std::vector<Fr>& a, const std::vector<Fr>& b);
    
    // Polynomial division (general)
    static std::vector<Fr> divide(const std::vector<Fr>& dividend, const std::vector<Fr>& divisor);
    
    // Optimized division by linear polynomial (x - z)
    static std::vector<Fr> divide_by_linear(const std::vector<Fr>& dividend, const Fr& z);
    
    // Lagrange interpolation for arbitrary points
    static std::vector<Fr> interpolate_lagrange(const std::vector<Fr>& x_vals, const std::vector<Fr>& y_vals);
    
    // NTT-based interpolation for roots of unity
    static std::vector<Fr> interpolate_ntt(const std::vector<Fr>& x_vals, const std::vector<Fr>& y_vals);
    
    // Generate random polynomial
    static std::vector<Fr> random(size_t degree);
    
    // Construct vanishing polynomial Z_H(x) = x^l - 1
    static std::vector<Fr> vanishing(size_t l);
    
    // Compute sum of evaluations on subgroup
    static Fr sum_on_subgroup(const std::vector<Fr>& poly, const Fr& omega, size_t l);
};