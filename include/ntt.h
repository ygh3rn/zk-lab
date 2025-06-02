#ifndef NTT_H
#define NTT_H

#include <vector>
#include <mcl/bn.hpp>

using namespace mcl::bn;

class NTT {
private:
    Fr omega;                   // primitive root of unity
    Fr omega_inv;               // inverse of omega
    size_t n;                   // transform size (must be power of 2)
    std::vector<Fr> powers;     // precomputed powers of omega
    std::vector<Fr> inv_powers; // precomputed powers of omega^(-1)
    
public:
    NTT(size_t size);
    
    // Forward NTT transform
    void forward_transform(std::vector<Fr>& a);
    
    // Inverse NTT transform
    void inverse_transform(std::vector<Fr>& a);
    
    // Utility functions
    Fr find_primitive_root(size_t n);
    void precompute_powers();
    bool is_primitive_root(const Fr& candidate, size_t n);
};

#endif