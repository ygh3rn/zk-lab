#ifndef KZG_H
#define KZG_H

#include <vector>
#include <mcl/bn.hpp>
#include "polynomial.h"

using namespace mcl::bn;

// Helper function for G2 commitments (used by both KZG and PIOP)
G2 commit_g2_internal(const Polynomial& poly, const std::vector<G2>& g2_powers);

class KZG {
private:
    std::vector<G1> g1_powers;  // [g₁, τg₁, τ²g₁, ..., τᵗg₁]
    std::vector<G2> g2_powers;  // [g₂, τg₂, τ²g₂, ..., τᵗg₂] - Extended for PIOP verification
    Fr tau;                     // secret parameter (for setup only)
    size_t max_degree;
    
public:
    KZG(size_t degree);
    
    // Setup phase
    void setup(size_t degree);
    
    // Commitment phase
    G1 commit(const Polynomial& poly);
    
    // Witness creation
    G1 create_witness(const Polynomial& poly, const Fr& point);
    
    // Verification
    bool verify_eval(const G1& commitment, const Fr& point, const Fr& value, const G1& witness);
    
    // Batch opening (optional)
    G1 create_batch_witness(const Polynomial& poly, const std::vector<Fr>& points);
    bool verify_batch_eval(const G1& commitment, const std::vector<Fr>& points, 
                          const std::vector<Fr>& values, const G1& witness);
    
    // Getters
    const std::vector<G1>& get_g1_powers() const { return g1_powers; }
    const std::vector<G2>& get_g2_powers() const { return g2_powers; }
};

#endif