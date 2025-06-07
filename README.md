# Cryptographic Protocol Implementation

A mathematically rigorous implementation of fundamental cryptographic protocols including polynomial commitments, Number Theoretic Transforms, and zero-knowledge Polynomial Interactive Oracle Proofs (PIOPs).

## Security & Mathematical Rigor

This implementation prioritizes cryptographic soundness and mathematical correctness. All protocols include:

- Formal security proofs in code comments
- Complete error handling with meaningful exceptions
- Input validation and bounds checking
- Production-ready verification with enhanced security checks
- Comprehensive testing against known attack vectors

## Features

### Core Cryptographic Primitives

#### 1. Number Theoretic Transform (NTT)
- Non-recursive implementation with optimized bit-reversal
- Mathematically verified primitive root finding for arbitrary field orders
- Time complexity: O(n log n) with minimal constant factors
- Applications: Fast polynomial multiplication, interpolation, evaluation

#### 2. KZG Polynomial Commitment Scheme
- Trusted setup with structured reference string generation
- Constant-size commitments (~32 bytes) for arbitrary degree polynomials
- Succinct evaluation proofs (~64 bytes) with pairing-based verification
- Knowledge soundness under the q-SDH assumption in bilinear groups

#### 3. Univariate ZeroTest PIOP
**Relation**: `R_UniZT = {(C ∈ 𝔾; f(X) ∈ 𝔽[X]) : C = PCS(f) ∧ ∀i ∈ [0,l), f(ωⁱ) = 0}`

- Mathematical foundation: Proves polynomial divisibility by vanishing polynomial Z_H(x) = x^l - 1
- Verification equation: e(C_f, g₂) = e(C_q, Z_H(τ)·g₂) using optimal ate pairing
- Security: Knowledge soundness assuming KZG binding and discrete log hardness

#### 4. Univariate SumCheck PIOP  
**Relation**: `R_UniSC = {(C ∈ 𝔾; f(X) ∈ 𝔽[X]) : C = PCS(f) ∧ ∑ᵢ₌₀ˡ⁻¹ f(ωⁱ) = 0}`

- Mathematical foundation: ∑f(ωⁱ) = l × [constant term of f(x) mod (x^l - 1)]
- Divisibility proof: Shows remainder r(x) = f(x) mod (x^l - 1) has r(x) = x·q(x)
- Pairing verification: e(C_f, g₂) = e(C_q, τ·g₂) for quotient witness

## Technical Specifications

### Cryptographic Parameters
- **Elliptic Curve**: BN_SNARK1 (254-bit prime order)
- **Security Level**: ~128 bits (suitable for production)
- **Pairing Type**: Optimal ate pairing over BN curves
- **Field Arithmetic**: Montgomery ladder with constant-time operations

### Performance Characteristics

| Operation | Time Complexity | Space Complexity | Proof Size |
|-----------|----------------|------------------|------------|
| NTT/INTT | O(n log n) | O(n) | N/A |
| KZG Setup | O(n)𝔾 | O(n) | O(n) SRS |
| KZG Commit | O(n)𝔾 | O(1) | 32 bytes |
| KZG Prove | O(n)𝔽 | O(n) | 64 bytes |
| KZG Verify | O(1)ₚ | O(1) | N/A |
| ZeroTest Prove | O(D)𝔾 + O(D)𝔽 | O(D) | 64 bytes |
| ZeroTest Verify | O(1)𝔾 + O(1)ₚ | O(1) | N/A |
| SumCheck Prove | O(D)𝔾 + O(D)𝔽 | O(D) | 64 bytes |
| SumCheck Verify | O(1)𝔾 + O(1)ₚ | O(1) | N/A |

*Legend: 𝔾 = group operation, 𝔽 = field operation, ₚ = pairing, D = polynomial degree*

## Architecture

```
include/
├── kzg.h              # KZG polynomial commitment interface
├── ntt.h              # Number theoretic transform operations  
├── polynomial.h       # Polynomial arithmetic and utilities
├── zerotest.h         # ZeroTest PIOP protocol
└── sumcheck.h         # SumCheck PIOP protocol

src/
├── kzg.cpp            # KZG implementation with security proofs
├── ntt.cpp            # Optimized NTT with primitive root finding
├── polynomial.cpp     # Polynomial operations (eval, multiply, divide)
├── zerotest.cpp       # ZeroTest prover/verifier with full verification
└── sumcheck.cpp       # SumCheck prover/verifier with mathematical rigor

tests/
└── test_suite.cpp     # Comprehensive test suite with attack vectors
```

## Getting Started

### Prerequisites

```bash
# Install MCL cryptographic library
git clone https://github.com/herumi/mcl.git
cd mcl && make -j$(nproc) && sudo make install

# System requirements
sudo apt update && sudo apt install -y cmake g++ libomp-dev
```

### Compilation

```bash
# Clone and build
git clone https://github.com/yourusername/cryptography-implementation.git
cd cryptography-implementation

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Run comprehensive test suite
./test_suite
```

### Basic Usage

```cpp
#include <mcl/bn.hpp>
#include "kzg.h"
#include "zerotest.h"
#include "sumcheck.h"

// Initialize curve parameters
mcl::initPairing(mcl::BN_SNARK1);

// Setup phase (trusted setup)
KZG::SetupParams params = KZG::Setup(1024);

// Create and commit to polynomial
std::vector<Fr> poly = {Fr(1), Fr(2), Fr(3)}; // 1 + 2x + 3x²
KZG::Commitment commit = KZG::Commit(poly, params);

// Generate and verify evaluation proof
Fr point = Fr(42);
KZG::Proof proof = KZG::CreateWitness(poly, point, params);
bool valid = KZG::VerifyEval(commit, point, proof, params);

// ZeroTest PIOP for vanishing polynomials
Fr omega = NTT::find_primitive_root(8);
std::vector<Fr> vanishing_poly = Polynomial::vanishing(8);
ZeroTestProof zero_proof = ZeroTest::prove(vanishing_poly, omega, 8, params);
bool zero_valid = ZeroTest::verify_with_full_checks(zero_proof, omega, 8, params);
```

## Testing & Verification

### Automated Test Suite

```bash
# Run all tests with performance benchmarks
./test_suite

# Expected output:
# NTT Implementation - correctness and timing
# Polynomial Operations - multiplication, division, interpolation  
# KZG Protocol - setup, commit, prove, verify
# ZeroTest PIOP - vanishing polynomial proofs
# SumCheck PIOP - sum-zero polynomial proofs
# Security Properties - completeness, soundness, binding
```

### Security Testing

The implementation includes tests for:
- **Completeness**: Honest provers always convince honest verifiers
- **Soundness**: Malicious provers cannot convince honest verifiers  
- **Binding**: Commitments cannot be opened to different values
- **Knowledge soundness**: Accepting provers can extract valid witnesses

## Mathematical Background

### Zero-Knowledge Polynomial IOPs

PIOPs combine the efficiency of Interactive Oracle Proofs with the strong security guarantees of zero-knowledge protocols. Key properties:

1. **Completeness**: If statement is true, honest prover convinces verifier
2. **Soundness**: If statement is false, no prover convinces verifier  
3. **Zero-Knowledge**: Verifier learns nothing beyond validity
4. **Succinctness**: Communication is sublinear in statement size

### KZG Commitments

Based on the mathematical insight that polynomial evaluation can be verified using bilinear pairings:

```
e(C - v·g₁, g₂) = e(π, τ·g₂ - z·g₂)
```

Where:
- `C = g₁^{f(τ)}` is polynomial commitment  
- `π = g₁^{(f(τ)-v)/(τ-z)}` is evaluation proof
- `v = f(z)` is claimed evaluation
- Security relies on q-Strong Diffie-Hellman assumption

## Security Considerations

### Cryptographic Assumptions
- **q-Strong Diffie-Hellman**: For KZG binding
- **Discrete Logarithm**: For commitment hiding  
- **Bilinear Diffie-Hellman**: For pairing security
- **Random Oracle Model**: For Fiat-Shamir transformation

### Implementation Security
- **Constant-time operations** prevent timing attacks
- **Input validation** prevents malformed data attacks
- **Memory safety** with RAII and bounds checking
- **Side-channel resistance** with uniform execution paths

### Trusted Setup
- KZG requires one-time trusted setup with toxic waste disposal
- Setup ceremony must be performed securely with multiple parties
- Public parameters can be reused across different applications
- Verification of setup correctness using mathematical properties

## References

1. **[KZG10]** Kate, A., Zaverucha, G. M., & Goldberg, I. (2010). *Constant-size commitments to polynomials and their applications*. ASIACRYPT 2010. Available at: https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf

2. **[Thaler22]** Thaler, J. (2022). *Proofs, Arguments, and Zero-Knowledge*. Available at: https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf

3. **MCL Library**: https://github.com/herumi/mcl - High-performance cryptographic library

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Contributions welcome! Please ensure:
- Mathematical correctness with formal verification
- Comprehensive testing including edge cases
- Security analysis for new protocols  
- Performance benchmarks for optimizations
- Code documentation with security proofs

---

*This implementation is for educational and research purposes. For production deployment, conduct thorough security audits and formal verification.*