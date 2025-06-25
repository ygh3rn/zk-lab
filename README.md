# Cryptographic Protocol Implementation

A mathematically rigorous implementation of fundamental cryptographic protocols including Number Theoretic Transforms, polynomial commitments, and zero-knowledge Polynomial Interactive Oracle Proofs (PIOPs).

## Number Theoretic Transform (NTT)

The Number Theoretic Transform is the finite field analogue of the Fast Fourier Transform, enabling efficient polynomial operations over finite fields.

### Mathematical Foundation

For prime *p* and primitive *n*-th root of unity *ω* where *n* | (*p*-1):

- **Forward NTT:** For polynomial *f(x)* = ∑ᵢ₌₀ⁿ⁻¹ *aᵢ xⁱ*:

<img src="https://latex.codecogs.com/svg.latex?\hat{f}_j%20=%20\sum_{i=0}^{n-1}%20a_i%20\omega^{ij}%20\bmod%20p%20\quad%20\text{for%20}%20j%20=%200,%20\ldots,%20n-1" />

- **Inverse NTT:** Reconstruction using *ω⁻¹*:

<img src="https://latex.codecogs.com/svg.latex?a_i%20=%20n^{-1}%20\sum_{j=0}^{n-1}%20\hat{f}_j%20(\omega^{-1})^{ij}%20\bmod%20p" />

### Implementation Details

Non-recursive implementation using bit-reversal permutation for optimal cache performance. Polynomial interpolation achieved through inverse NTT of evaluation points.

## Polynomial Multiplication via NTT

Efficient multiplication using convolution theorem in finite fields.

### Algorithm

For polynomials *f(x)*, *g(x)* of degree less than *n*:

1. Pad to size 2*n* and compute NTT(*f*), NTT(*g*)
2. Pointwise multiplication: NTT(*h*)[*i*] = NTT(*f*)[*i*] · NTT(*g*)[*i*]
3. Inverse transform: *h(x)* = INTT(NTT(*h*))

### Complexity Analysis

- **Time complexity:** O(*n* log *n*) versus naive O(*n*²)
- **Space complexity:** O(*n*) field elements

## KZG Polynomial Commitment Scheme

Implementation of the Kate-Zaverucha-Goldberg commitment scheme for univariate polynomials.

### Interface

- **Setup**(1^λ, *d*) → SRS: Generate (*g*, *g^τ*, *g^τ²*, ..., *g^τᵈ*)
- **Commit**(*f(x)*) → *c*: Compute *c* = *g^f(τ)*
- **CreateWitness**(*f(x)*, *z*) → *π*: Generate *π* = *g^w(τ)* where *w(x)* = (*f(x)* - *f(z)*)/(*x* - *z*)
- **VerifyEval**(*c*, *z*, *v*, *π*) → {0,1}: Check *e(c · g^(-v), g)* = *e(π, g^τ · g^(-z))*

### Security Properties

- **Binding:** Under *d*-Strong Diffie-Hellman assumption
- **Succinctness:** Constant proof size independent of polynomial degree

## Univariate ZeroTest PIOP

A PIOP proving that polynomial *f(x)* evaluates to zero everywhere on subgroup 𝐇\_ℓ.

### Relation

<img src="https://latex.codecogs.com/svg.latex?\mathcal{R}_{\text{UniZT}}%20=%20\{(C%20\in%20\mathbb{G};%20f(X)%20\in%20\mathbb{F}_p[X])%20:%20C%20=%20\text{PCS}(f)%20\land%20\forall%20i%20\in%20[0,\ell),%20f(\omega_\ell^i)%20=%200\}" />

### Protocol Description

**Key insight:** *f(x)* = 0 on 𝐇\_ℓ ⟺ Z₍𝐇\_ℓ₎(*x*) | *f(x)* where Z₍𝐇\_ℓ₎(*x*) = ∏ᵢ₌₀^(ℓ-1)(*x* - ωₗⁱ).

1. **Prover:** Compute quotient *q(x)* = *f(x)*/Z₍𝐇\_ℓ₎(*x*) and send *Cᵩ* = PCS(*q*)
2. **Verifier:** Send random challenge *r* ∈ 𝔽\_p
3. **Prover:** Open *f(r)* and *q(r)* with respective proofs
4. **Verifier:** Check *f(r)* = *q(r)* · Z₍𝐇\_ℓ₎(*r*) and verify opening proofs

### Complexity Analysis

- **Prover time:** O(*D*)𝔾 + O(*D*)𝔽
- **Verifier time:** O(1)𝔾 + O(1)𝔽
- **Proof size:** O(1)

## Univariate SumCheck PIOP

A PIOP proving that polynomial evaluations sum to zero over subgroup 𝐇\_ℓ.

### Relation

<img src="https://latex.codecogs.com/svg.latex?\mathcal{R}_{\text{UniSC}}%20=%20\left\{(C%20\in%20\mathbb{G};%20f(X)%20\in%20\mathbb{F}_p[X])%20:%20C%20=%20\text{PCS}(f)%20\land%20\sum_{i=0}^{\ell-1}%20f(\omega_\ell^i)%20=%200\right\}" />

### Protocol Description

**Key insight:** For polynomial *g(x)* of degree < ℓ: ∑\_(*a* ∈ 𝐇\_ℓ) *g(a)* = *g(0)* · ℓ.

Therefore, ∑\_(*a* ∈ 𝐇\_ℓ) *f(a)* = 0 ⟺ ∃ h\*(*x*), *t(x)* such that:

<img src="https://latex.codecogs.com/svg.latex?f(x)%20=%20h^*(x)%20\cdot%20Z_{\mathbb{H}_\ell}(x)%20+%20x%20\cdot%20t(x)" />

where deg(*t*) < ℓ - 1.

1. **Prover:** Compute h\*(*x*) and *t(x)*, send C\_{h\*} = PCS(h\*) and *Cₜ* = PCS(*t*)
2. **Verifier:** Send random challenge *r* ∈ 𝔽\_p
3. **Prover:** Open *f(r)*, h\*(*r*), and *t(r)* with respective proofs
4. **Verifier:** Check *f(r)* = h\*(*r*) · Z₍𝐇\_ℓ₎(*r*) + *r* · *t(r)* and verify opening proofs

### Complexity Analysis

- **Prover time:** O(*D*)𝔾 + O(*D*)𝔽
- **Verifier time:** O(1)𝔾 + O(1)𝔽
- **Proof size:** O(1)

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
git clone <your-repo-url>
cd cryptography-implementation

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)

# Run comprehensive test suite
./test_suite
```

### Project Structure

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

## Implementation Notes

- **Library:** Used mcl library with BN\_SNARK1 curve for cryptographic operations
- **Testing:** Generated random polynomials and verified against mathematical properties
- **Security:** All protocols include formal verification and security proofs
- **Performance:** Optimized implementations with complexity analysis

## Usage Example

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
```

## References

1. **[KZG10]** Kate, A., Zaverucha, G. M., & Goldberg, I. (2010). *Constant-size commitments to polynomials and their applications*. ASIACRYPT 2010.

2. **[Thaler22]** Thaler, J. (2022). *Proofs, Arguments, and Zero-Knowledge*.

3. **MCL Library**: https://github.com/herumi/mcl - High-performance cryptographic library

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

*This implementation is for educational and research purposes. Mathematical rigor and security are prioritized throughout.*