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
    
    // Basic operations
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;  // Using NTT
    
    // Evaluation
    Fr evaluate(const Fr& x) const;
    
    // Polynomial interpolation
    static Polynomial lagrange_interpolation(const std::vector<std::pair<Fr, Fr>>& points);
    
    // Getters/Setters
    const std::vector<Fr>& get_coefficients() const { return coefficients; }
    void set_coefficients(const std::vector<Fr>& coeffs) { coefficients = coeffs; }
    size_t degree() const;
    
    // Utility functions
    void resize(size_t new_size);
    void print() const;
};

#endif