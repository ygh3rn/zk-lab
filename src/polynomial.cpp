#include "polynomial.h"
#include <algorithm>
#include <iostream>
#include <cassert>
#include <stdexcept>

Polynomial::Polynomial(const std::vector<Fr>& coeffs) : coefficients(coeffs) {
    remove_leading_zeros();
    if (coefficients.empty()) {
        coefficients.push_back(Fr(0));
    }
}

void Polynomial::remove_leading_zeros() {
    while (coefficients.size() > 1 && coefficients.back().isZero()) {
        coefficients.pop_back();
    }
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    size_t max_size = std::max(coefficients.size(), other.coefficients.size());
    std::vector<Fr> result(max_size);
    
    for (size_t i = 0; i < max_size; ++i) {
        Fr a = (i < coefficients.size()) ? coefficients[i] : Fr(0);
        Fr b = (i < other.coefficients.size()) ? other.coefficients[i] : Fr(0);
        result[i] = a + b;
    }
    
    return Polynomial(result);
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    size_t max_size = std::max(coefficients.size(), other.coefficients.size());
    std::vector<Fr> result(max_size);
    
    for (size_t i = 0; i < max_size; ++i) {
        Fr a = (i < coefficients.size()) ? coefficients[i] : Fr(0);
        Fr b = (i < other.coefficients.size()) ? other.coefficients[i] : Fr(0);
        result[i] = a - b;
    }
    
    return Polynomial(result);
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    if (is_zero() || other.is_zero()) {
        return Polynomial({Fr(0)});
    }
    
    // Determine the size of the result
    size_t result_degree = degree() + other.degree();
    size_t n = 1;
    while (n <= result_degree) {
        n <<= 1;
    }
    n <<= 1; // Double for padding to avoid wrap-around
    
    // Create NTT transformer
    NTT ntt(n);
    
    // Prepare padded inputs
    std::vector<Fr> a(n, Fr(0));
    std::vector<Fr> b(n, Fr(0));
    
    std::copy(coefficients.begin(), coefficients.end(), a.begin());
    std::copy(other.coefficients.begin(), other.coefficients.end(), b.begin());
    
    // Forward NTT
    ntt.forward_transform(a);
    ntt.forward_transform(b);
    
    // Pointwise multiplication
    for (size_t i = 0; i < n; ++i) {
        a[i] *= b[i];
    }
    
    // Inverse NTT
    ntt.inverse_transform(a);
    
    // Extract result (remove padding)
    std::vector<Fr> result(result_degree + 1);
    std::copy(a.begin(), a.begin() + result_degree + 1, result.begin());
    
    return Polynomial(result);
}

// MAIN DIVISION OPERATOR - chooses optimal algorithm automatically
Polynomial Polynomial::operator/(const Polynomial& divisor) const {
    if (divisor.is_zero()) {
        throw std::runtime_error("Division by zero polynomial");
    }
    
    // Optimization: For linear divisors (x - a), use synthetic division - O(d)
    if (divisor.coefficients.size() == 2 && divisor.coefficients[1] == Fr(1)) {
        Fr root = -divisor.coefficients[0]; // divisor is (x - root)
        return divide_by_linear_synthetic(root);
    }
    
    // Fall back to general long division - O(d²)
    return divide_general(divisor);
}

// O(d) synthetic division for (x - root)
Polynomial Polynomial::divide_by_linear_synthetic(const Fr& root) const {
    if (coefficients.size() <= 1) {
        return Polynomial({Fr(0)});
    }
    
    std::vector<Fr> quotient_coeffs(coefficients.size() - 1);
    Fr accumulator = coefficients.back();
    
    for (int i = static_cast<int>(coefficients.size()) - 2; i >= 0; --i) {
        quotient_coeffs[i] = accumulator;
        accumulator = accumulator * root + coefficients[i];
    }
    
    return Polynomial(quotient_coeffs);
}

// General polynomial division - O(d²)
Polynomial Polynomial::divide_general(const Polynomial& divisor) const {
    const auto& dividend_coeffs = coefficients;
    const auto& divisor_coeffs = divisor.coefficients;
    
    if (dividend_coeffs.size() < divisor_coeffs.size()) {
        return Polynomial({Fr(0)}); // quotient is zero
    }
    
    std::vector<Fr> quotient_coeffs(dividend_coeffs.size() - divisor_coeffs.size() + 1, Fr(0));
    std::vector<Fr> remainder_coeffs = dividend_coeffs;
    
    Fr leading_coeff_inv;
    Fr::inv(leading_coeff_inv, divisor_coeffs.back());
    
    for (int i = static_cast<int>(quotient_coeffs.size()) - 1; i >= 0; --i) {
        if (remainder_coeffs.size() >= divisor_coeffs.size()) {
            quotient_coeffs[i] = remainder_coeffs.back() * leading_coeff_inv;
            
            for (size_t j = 0; j < divisor_coeffs.size(); ++j) {
                size_t remainder_idx = remainder_coeffs.size() - divisor_coeffs.size() + j;
                remainder_coeffs[remainder_idx] -= quotient_coeffs[i] * divisor_coeffs[j];
            }
            
            remainder_coeffs.pop_back();
        }
    }
    
    return Polynomial(quotient_coeffs);
}

// Create vanishing polynomial Z_H(x) = x^n - 1
Polynomial Polynomial::vanishing_polynomial(size_t subgroup_size) {
    std::vector<Fr> coeffs(subgroup_size + 1, Fr(0));
    coeffs[0] = Fr(-1);           // -1
    coeffs[subgroup_size] = Fr(1); // x^n
    return Polynomial(coeffs);
}

// Create polynomial from roots ∏(x - root_i)
Polynomial Polynomial::from_roots(const std::vector<Fr>& roots) {
    Polynomial result({Fr(1)}); // Start with polynomial 1
    
    for (const Fr& root : roots) {
        Polynomial linear_factor({-root, Fr(1)}); // (x - root)
        result = result * linear_factor;
    }
    
    return result;
}

std::pair<Polynomial, Polynomial> Polynomial::divmod(const Polynomial& divisor) const {
    if (divisor.is_zero()) {
        throw std::runtime_error("Division by zero polynomial");
    }
    
    // Check for special cases first (optimizations)
    if (divisor.degree() == 1 && divisor.coefficients[1] == Fr(1)) {
        // Linear divisor (x - a) - use O(d) synthetic division
        Fr root = -divisor.coefficients[0];
        Polynomial quotient = divide_by_linear(root);
        Fr remainder = evaluate(root);
        return {quotient, Polynomial({remainder})};
    }
    
    // Fall back to general long division
    return long_division(divisor);
}

// O(d) synthetic division for (x - root)
Polynomial Polynomial::divide_by_linear(const Fr& root) const {
    if (coefficients.size() <= 1) {
        return Polynomial({Fr(0)});
    }
    
    std::vector<Fr> quotient_coeffs(coefficients.size() - 1);
    Fr accumulator = coefficients.back();
    
    for (int i = static_cast<int>(coefficients.size()) - 2; i >= 0; --i) {
        quotient_coeffs[i] = accumulator;
        accumulator = accumulator * root + coefficients[i];
    }
    // Accumulator is now the remainder (should equal evaluate(root))
    
    return Polynomial(quotient_coeffs);
}

// Division by vanishing polynomial x^n - 1
Polynomial Polynomial::divide_by_vanishing(size_t subgroup_size) const {
    if (coefficients.size() <= subgroup_size) {
        return Polynomial({Fr(0)});
    }
    
    // For x^n - 1, the quotient has a simple structure
    std::vector<Fr> quotient_coeffs(coefficients.size() - subgroup_size);
    
    for (size_t i = 0; i < quotient_coeffs.size(); ++i) {
        quotient_coeffs[i] = coefficients[i + subgroup_size];
        if (i < coefficients.size() - subgroup_size) {
            quotient_coeffs[i] += coefficients[i]; // Add lower degree terms
        }
    }
    
    return Polynomial(quotient_coeffs);
}

// General polynomial long division - O(d²) but works for any divisor
std::pair<Polynomial, Polynomial> Polynomial::long_division(const Polynomial& divisor) const {
    const auto& dividend_coeffs = coefficients;
    const auto& divisor_coeffs = divisor.coefficients;
    
    if (dividend_coeffs.size() < divisor_coeffs.size()) {
        return {Polynomial({Fr(0)}), *this}; // quotient=0, remainder=dividend
    }
    
    std::vector<Fr> quotient_coeffs(dividend_coeffs.size() - divisor_coeffs.size() + 1, Fr(0));
    std::vector<Fr> remainder_coeffs = dividend_coeffs;
    
    Fr leading_coeff_inv;
    Fr::inv(leading_coeff_inv, divisor_coeffs.back());
    
    for (int i = static_cast<int>(quotient_coeffs.size()) - 1; i >= 0; --i) {
        if (remainder_coeffs.size() >= divisor_coeffs.size()) {
            // Compute quotient coefficient
            quotient_coeffs[i] = remainder_coeffs.back() * leading_coeff_inv;
            
            // Subtract divisor * quotient_coeff from remainder
            for (size_t j = 0; j < divisor_coeffs.size(); ++j) {
                size_t remainder_idx = remainder_coeffs.size() - divisor_coeffs.size() + j;
                remainder_coeffs[remainder_idx] -= quotient_coeffs[i] * divisor_coeffs[j];
            }
            
            // Remove leading zero
            remainder_coeffs.pop_back();
        }
    }
    
    return {Polynomial(quotient_coeffs), Polynomial(remainder_coeffs)};
}

Fr Polynomial::evaluate(const Fr& x) const {
    if (coefficients.empty()) {
        return Fr(0);
    }
    
    // Horner's method
    Fr result = coefficients.back();
    for (int i = static_cast<int>(coefficients.size()) - 2; i >= 0; --i) {
        result = result * x + coefficients[i];
    }
    
    return result;
}

Polynomial Polynomial::lagrange_interpolation(const std::vector<std::pair<Fr, Fr>>& points) {
    size_t n = points.size();
    if (n == 0) {
        return Polynomial({Fr(0)});
    }
    
    std::vector<Fr> result_coeffs(n, Fr(0));
    
    for (size_t i = 0; i < n; ++i) {
        // Compute the i-th Lagrange basis polynomial
        std::vector<Fr> basis_coeffs = {Fr(1)};
        
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                Fr xi = points[i].first;
                Fr xj = points[j].first;
                Fr denominator = xi - xj;
                Fr denominator_inv;
                Fr::inv(denominator_inv, denominator);
                
                // Multiply basis polynomial by (x - xj) / (xi - xj)
                std::vector<Fr> new_coeffs(basis_coeffs.size() + 1, Fr(0));
                
                for (size_t k = 0; k < basis_coeffs.size(); ++k) {
                    new_coeffs[k] += basis_coeffs[k] * (-xj) * denominator_inv;
                    new_coeffs[k + 1] += basis_coeffs[k] * denominator_inv;
                }
                
                basis_coeffs = new_coeffs;
            }
        }
        
        // Scale by yi and add to result
        Fr yi = points[i].second;
        for (size_t k = 0; k < basis_coeffs.size() && k < result_coeffs.size(); ++k) {
            result_coeffs[k] += basis_coeffs[k] * yi;
        }
    }
    
    return Polynomial(result_coeffs);
}

size_t Polynomial::degree() const {
    if (coefficients.size() <= 1) {
        return coefficients.empty() || coefficients[0].isZero() ? 0 : 0;
    }
    return coefficients.size() - 1;
}

bool Polynomial::is_zero() const {
    return coefficients.size() == 1 && coefficients[0].isZero();
}

void Polynomial::resize(size_t new_size) {
    coefficients.resize(new_size, Fr(0));
}

void Polynomial::print() const {
    std::cout << "Polynomial: ";
    bool first = true;
    
    for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; --i) {
        if (!coefficients[i].isZero()) {
            if (!first && !coefficients[i].isZero()) {
                std::cout << " + ";
            }
            first = false;
            
            if (i == 0) {
                std::cout << coefficients[i];
            } else if (i == 1) {
                if (coefficients[i] == Fr(1)) {
                    std::cout << "x";
                } else {
                    std::cout << coefficients[i] << "*x";
                }
            } else {
                if (coefficients[i] == Fr(1)) {
                    std::cout << "x^" << i;
                } else {
                    std::cout << coefficients[i] << "*x^" << i;
                }
            }
        }
    }
    
    if (first) {
        std::cout << "0";
    }
    
    std::cout << std::endl;
}