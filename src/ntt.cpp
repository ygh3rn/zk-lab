#include "ntt.h"
#include <algorithm>
#include <cassert>
#include <iostream>

NTT::NTT(size_t size) : n(size) {
    // Ensure n is a power of 2
    assert((n & (n - 1)) == 0 && "Size must be a power of 2");
    
    // Find primitive n-th root of unity
    omega = find_primitive_root(n);
    
    // Compute inverse of omega
    Fr::inv(omega_inv, omega);
    
    // Precompute powers
    precompute_powers();
}

Fr NTT::find_primitive_root(size_t n) {
    // For BN curves, we know that primitive roots exist 
    Fr candidate;
    
    for (int i = 2; i < 100; ++i) {
        candidate = Fr(i);
        if (is_primitive_root(candidate, n)) {
            return candidate;
        }
    }
    
    // If no small primitive root found, try a different approach
    // Use the fact that g^((p-1)/n) gives us an n-th root of unity
    // where g is a generator of the multiplicative group
    
    candidate = Fr(5); // Known to be a generator in many fields
    
    // Compute candidate^((field_order-1)/n)
    // We'll approximate this for the BN curve
    Fr exponent;
    exponent = Fr(-1);
    exponent /= Fr(n);
    
    Fr::pow(candidate, candidate, exponent);
    
    return candidate;
}

bool NTT::is_primitive_root(const Fr& candidate, size_t n) {
    Fr temp = candidate;
    
    // Check if candidate^n ≡ 1
    Fr::pow(temp, candidate, Fr(n));
    if (!temp.isOne()) {
        return false;
    }
    
    // Check that candidate^k ≢ 1 for proper divisors k of n
    for (size_t k = 1; k < n; ++k) {
        if (n % k == 0) {
            Fr::pow(temp, candidate, Fr(k));
            if (temp.isOne()) {
                return false;
            }
        }
    }
    
    return true;
}

void NTT::precompute_powers() {
    powers.resize(n);
    inv_powers.resize(n);
    
    powers[0] = Fr(1);
    inv_powers[0] = Fr(1);
    
    for (size_t i = 1; i < n; ++i) {
        powers[i] = powers[i-1] * omega;
        inv_powers[i] = inv_powers[i-1] * omega_inv;
    }
}

void NTT::forward_transform(std::vector<Fr>& a) {
    assert(a.size() == n && "Input size must match NTT size");
    
    // Bit-reversal permutation
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
    
    // Cooley-Tukey NTT
    for (size_t length = 2; length <= n; length <<= 1) {
        Fr wlen = powers[n / length];
        
        for (size_t i = 0; i < n; i += length) {
            Fr w(1);
            
            for (size_t j = 0; j < length / 2; ++j) {
                Fr u = a[i + j];
                Fr v = a[i + j + length / 2] * w;
                
                a[i + j] = u + v;
                a[i + j + length / 2] = u - v;
                
                w *= wlen;
            }
        }
    }
}

void NTT::inverse_transform(std::vector<Fr>& a) {
    assert(a.size() == n && "Input size must match NTT size");
    
    // Bit-reversal permutation
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
    
    // Cooley-Tukey INTT (using inverse roots)
    for (size_t length = 2; length <= n; length <<= 1) {
        Fr wlen = inv_powers[n / length];
        
        for (size_t i = 0; i < n; i += length) {
            Fr w(1);
            
            for (size_t j = 0; j < length / 2; ++j) {
                Fr u = a[i + j];
                Fr v = a[i + j + length / 2] * w;
                
                a[i + j] = u + v;
                a[i + j + length / 2] = u - v;
                
                w *= wlen;
            }
        }
    }
    
    // Divide by n
    Fr n_inv;
    Fr::inv(n_inv, Fr(n));
    
    for (size_t i = 0; i < n; ++i) {
        a[i] *= n_inv;
    }
}