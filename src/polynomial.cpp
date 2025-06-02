#include "polynomial.h"
#include <algorithm>
#include <iostream>
#include <cassert>

Polynomial::Polynomial(const std::vector<Fr>& coeffs) : coefficients(coeffs) {
    // Remove leading zeros
    while (coefficients.size() > 1 && coefficients.back().isZero()) {
        coefficients.pop_back();
    }
    
    if (coefficients.empty()) {
        coefficients.push_back(Fr(0));
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
    if (coefficients.size() == 1 && coefficients[0].isZero()) {
        return Polynomial({Fr(0)});
    }
    if (other.coefficients.size() == 1 && other.coefficients[0].isZero()) {
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