#pragma once

#include <mcl/bn.hpp>
#include <vector>

using namespace mcl;
using namespace std;

class NTT {
public:
    // Non-recursive NTT with bit-reversal
    static vector<Fr> transform(const vector<Fr>& a, const Fr& root, size_t n);
    
    // Inverse NTT
    static vector<Fr> inverse_transform(const vector<Fr>& a, const Fr& root, size_t n);
    
    // Find primitive nth root of unity
    static Fr find_primitive_root(size_t n);
    
private:
    static Fr mod_inverse(const Fr& a);
};