#include "ntt.h"
#include <stdexcept>
#include <algorithm>

vector<Fr> NTT::transform(const vector<Fr>& a, const Fr& root, size_t n) {
    if (n == 0 || (n & (n - 1)) != 0) {
        throw invalid_argument("NTT size must be power of 2");
    }
    
    vector<Fr> result = a;
    result.resize(n, Fr(0));
    
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            swap(result[i], result[j]);
        }
    }
    
    for (size_t len = 2; len <= n; len <<= 1) {
        Fr wlen;
        Fr exponent = Fr(n / len);
        Fr::pow(wlen, root, exponent);
        
        for (size_t i = 0; i < n; i += len) {
            Fr w = Fr(1);
            
            for (size_t j = 0; j < len / 2; j++) {
                Fr u = result[i + j];
                Fr v;
                Fr::mul(v, result[i + j + len / 2], w);
                
                Fr::add(result[i + j], u, v);
                Fr::sub(result[i + j + len / 2], u, v);
                Fr::mul(w, w, wlen);
            }
        }
    }
    
    return result;
}

vector<Fr> NTT::inverse_transform(const vector<Fr>& a, const Fr& root, size_t n) {
    Fr inv_root = mod_inverse(root);
    vector<Fr> result = transform(a, inv_root, n);
    
    Fr inv_n = mod_inverse(Fr(n));
    for (auto& x : result) {
        Fr::mul(x, x, inv_n);
    }
    
    return result;
}

Fr NTT::find_primitive_root(size_t n) {
    if (n == 0 || (n & (n - 1)) != 0) {
        throw invalid_argument("n must be power of 2");
    }
    
    Fr field_order_minus_one = Fr(-1);
    Fr exponent;
    Fr::div(exponent, field_order_minus_one, Fr(n));
    
    vector<int> candidates = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    
    for (int base : candidates) {
        Fr candidate = Fr(base);
        Fr root;
        Fr::pow(root, candidate, exponent);
        
        Fr test_n;
        Fr::pow(test_n, root, Fr(n));
        
        if (test_n == Fr(1)) {
            bool is_primitive = true;
            
            for (size_t divisor = 2; divisor < n; divisor++) {
                if (n % divisor == 0) {
                    Fr test_order;
                    Fr::pow(test_order, root, Fr(divisor));
                    if (test_order == Fr(1)) {
                        is_primitive = false;
                        break;
                    }
                }
            }
                    
            if (is_primitive) {
                return root;
            }
        }
    }
    
    throw runtime_error("Failed to find primitive nth root of unity");
}

Fr NTT::mod_inverse(const Fr& a) {
    Fr result;
    Fr::inv(result, a);
    return result;
}