#include "polynomial.h"
#include "ntt.h"
#include <stdexcept>
#include <cassert>

Fr Polynomial::evaluate(const std::vector<Fr>& coefficients, const Fr& x) {
    if (coefficients.empty()) {
        return Fr(0);
    }
    
    Fr result = coefficients.back();
    for (int i = coefficients.size() - 2; i >= 0; i--) {
        Fr::mul(result, result, x);
        Fr::add(result, result, coefficients[i]);
    }
    
    return result;
}

std::vector<Fr> Polynomial::multiply(const std::vector<Fr>& a, const std::vector<Fr>& b) {
    if (a.empty() || b.empty()) {
        return {};
    }
    
    size_t result_size = a.size() + b.size() - 1;
    size_t n = 1;
    while (n < result_size) n <<= 1;
    
    Fr root = NTT::find_primitive_root(n);
    
    std::vector<Fr> fa = NTT::transform(a, root, n);
    std::vector<Fr> fb = NTT::transform(b, root, n);
    
    for (size_t i = 0; i < n; i++) {
        Fr::mul(fa[i], fa[i], fb[i]);
    }
    
    std::vector<Fr> result = NTT::inverse_transform(fa, root, n);
    result.resize(result_size);
    
    return result;
}

std::vector<Fr> Polynomial::divide_by_linear(const std::vector<Fr>& dividend, const Fr& z) {
    if (dividend.empty()) {
        return {};
    }
    
    if (dividend.size() == 1) {
        return {};
    }
    
    std::vector<Fr> quotient(dividend.size() - 1);
    quotient[quotient.size() - 1] = dividend.back();
    
    for (int i = quotient.size() - 2; i >= 0; i--) {
        Fr temp;
        Fr::mul(temp, quotient[i + 1], z);
        Fr::add(quotient[i], dividend[i + 1], temp);
    }
    
    return quotient;
}

std::vector<Fr> Polynomial::divide(const std::vector<Fr>& dividend, const std::vector<Fr>& divisor) {
    if (divisor.empty()) {
        throw std::invalid_argument("Divisor cannot be empty");
    }
    
    std::vector<Fr> clean_divisor = divisor;
    while (!clean_divisor.empty() && clean_divisor.back().isZero()) {
        clean_divisor.pop_back();
    }
    
    if (clean_divisor.empty() || clean_divisor.back().isZero()) {
        throw std::invalid_argument("Divisor cannot be zero polynomial");
    }
    
    std::vector<Fr> clean_dividend = dividend;
    while (!clean_dividend.empty() && clean_dividend.back().isZero()) {
        clean_dividend.pop_back();
    }
    
    if (clean_dividend.empty()) {
        return {};
    }
    
    if (clean_dividend.size() < clean_divisor.size()) {
        return {};
    }
    
    size_t quotient_degree = clean_dividend.size() - clean_divisor.size();
    std::vector<Fr> quotient(quotient_degree + 1, Fr(0));
    std::vector<Fr> remainder = clean_dividend;
    
    Fr leading_coeff_inv;
    Fr::inv(leading_coeff_inv, clean_divisor.back());
    
    for (int i = quotient_degree; i >= 0; i--) {
        if (remainder.size() >= clean_divisor.size()) {
            Fr::mul(quotient[i], remainder.back(), leading_coeff_inv);
            
            for (size_t j = 0; j < clean_divisor.size(); j++) {
                if (remainder.size() >= clean_divisor.size() - j) {
                    size_t remainder_idx = remainder.size() - clean_divisor.size() + j;
                    Fr term;
                    Fr::mul(term, quotient[i], clean_divisor[j]);
                    Fr::sub(remainder[remainder_idx], remainder[remainder_idx], term);
                }
            }
            
            if (!remainder.empty() && remainder.back().isZero()) {
                remainder.pop_back();
            }
        }
    }
    
    return quotient;
}

std::vector<Fr> Polynomial::interpolate_lagrange(const std::vector<Fr>& x_vals, const std::vector<Fr>& y_vals) {
    assert(x_vals.size() == y_vals.size());
    size_t n = x_vals.size();
    
    if (n == 0) return {};
    if (n == 1) return {y_vals[0]};
    
    std::vector<Fr> result(n, Fr(0));
    
    for (size_t i = 0; i < n; i++) {
        std::vector<Fr> basis = {Fr(1)};
        
        for (size_t j = 0; j < n; j++) {
            if (i != j) {
                Fr denominator;
                Fr::sub(denominator, x_vals[i], x_vals[j]);
                Fr inv_denom;
                Fr::inv(inv_denom, denominator);
                
                std::vector<Fr> linear = {Fr(0), Fr(1)};
                Fr::sub(linear[0], linear[0], x_vals[j]);
                
                basis = multiply(basis, linear);
                
                for (auto& coeff : basis) {
                    Fr::mul(coeff, coeff, inv_denom);
                }
            }
        }
        
        for (size_t k = 0; k < basis.size() && k < result.size(); k++) {
            Fr term;
            Fr::mul(term, basis[k], y_vals[i]);
            Fr::add(result[k], result[k], term);
        }
    }
    
    return result;
}

std::vector<Fr> Polynomial::interpolate_ntt(const std::vector<Fr>& x_vals, const std::vector<Fr>& y_vals) {
    size_t n = x_vals.size();
    if (n == 0) return {};
    
    if (n > 1) {
        Fr omega = NTT::find_primitive_root(n);
        Fr omega_power = Fr(1);
        
        bool is_subgroup = true;
        for (size_t i = 0; i < n; i++) {
            if (!(x_vals[i] == omega_power)) {
                is_subgroup = false;
                break;
            }
            Fr::mul(omega_power, omega_power, omega);
        }
        
        if (is_subgroup) {
            return NTT::inverse_transform(y_vals, omega, n);
        }
    }
    
    return interpolate_lagrange(x_vals, y_vals);
}

std::vector<Fr> Polynomial::random(size_t degree) {
    std::vector<Fr> poly(degree + 1);
    for (auto& coeff : poly) {
        coeff.setByCSPRNG();
    }
    return poly;
}

std::vector<Fr> Polynomial::vanishing(size_t l) {
    std::vector<Fr> result(l + 1, Fr(0));
    result[0] = Fr(-1);
    result[l] = Fr(1);
    return result;
}

Fr Polynomial::sum_on_subgroup(const std::vector<Fr>& poly, const Fr& omega, size_t l) {
    Fr sum = Fr(0);
    Fr omega_power = Fr(1);
    
    for (size_t i = 0; i < l; i++) {
        Fr eval = evaluate(poly, omega_power);
        Fr::add(sum, sum, eval);
        Fr::mul(omega_power, omega_power, omega);
    }
    
    return sum;
}