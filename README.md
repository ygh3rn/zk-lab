# Cryptographic Protocol Implementation

A mathematically rigorous implementation of fundamental cryptographic protocols including Number Theoretic Transforms, polynomial commitments, and zero-knowledge Polynomial Interactive Oracle Proofs (PIOPs).

## Number Theoretic Transform (NTT)

The Number Theoretic Transform is the finite field analogue of the Fast Fourier Transform, enabling efficient polynomial operations over finite fields.

### Mathematical Foundation

For prime _p_ and primitive _n_-th root of unity _ω_ where _n_ | (_p_-1):

- **Forward NTT:** For polynomial _f(x)_ = ∑ᵢ₌₀ⁿ⁻¹ _aᵢ xⁱ_:

<img src="https://latex.codecogs.com/svg.latex?\hat{f}_j%20=%20\sum_{i=0}^{n-1}%20a_i%20\omega^{ij}%20\bmod%20p%20\quad%20\text{for%20}%20j%20=%200,%20\ldots,%20n-1" />

- **Inverse NTT:** Reconstruction using _ω⁻¹_:

<img src="https://latex.codecogs.com/svg.latex?a_i%20=%20n^{-1}%20\sum_{j=0}^{n-1}%20\hat{f}_j%20(\omega^{-1})^{ij}%20\bmod%20p" />

### Implementation Details

Non-recursive implementation using bit-reversal permutation for optimal cache performance. Polynomial interpolation achieved through inverse NTT of evaluation points.

## Polynomial Multiplication via NTT

Efficient multiplication using convolution theorem in finite fields.

### Algorithm

For polynomials _f(x)_, _g(x)_ of degree less than _n_:

1. Pad to size 2_n_ and compute NTT(_f_), NTT(_g_)
2. Pointwise multiplication: NTT(_h_)[_i_] = NTT(_f_)[_i_] · NTT(_g_)[_i_]
3. Inverse transform: _h(x)_ = INTT(NTT(_h_))

### Complexity Analysis

- **Time complexity:** _O(n log n)_ versus naive _O(n²)_
- **Space complexity:** _O(n)_ field elements

## KZG Polynomial Commitment Scheme

Implementation of the Kate-Zaverucha-Goldberg commitment scheme for univariate polynomials.

### Interface

- **Setup**(1^λ, _d_) → SRS: Generate (_g_, _g^τ_, _g^τ²_, ..., _g^τᵈ_)
- **Commit**(_f(x)_) → _c_: Compute _c_ = _g^f(τ)_
- **CreateWitness**(_f(x)_, _z_) → _π_: Generate _π_ = _g^w(τ)_ where _w(x)_ = (_f(x)_ - _f(z)_)/(_x_ - _z_)
- **VerifyEval**(_c_, _z_, _v_, _π_) → {0,1}: Check _e(c · g^(-v), g)_ = _e(π, g^τ · g^(-z))_

### Security Properties

- **Binding:** Under _d_-Strong Diffie-Hellman assumption
- **Succinctness:** Constant proof size independent of polynomial degree

## Univariate ZeroTest PIOP

A PIOP proving that polynomial _f(x)_ evaluates to zero everywhere on subgroup 𝐇_ℓ.

### Relation

<img src="https://latex.codecogs.com/svg.latex?\mathcal{R}_{\text{UniZT}}%20=%20\{(C%20\in%20\mathbb{G};%20f(X)%20\in%20\mathbb{F}_p[X])%20:%20C%20=%20\text{PCS}(f)%20\land%20\forall%20i%20\in%20[0,\ell),%20f(\omega_\ell^i)%20=%200\}" />

### Protocol Description

**Key insight:** _f(x)_ = 0 on 𝐇_ℓ ⟺ Z₍𝐇_ℓ₎(_x_) | _f(x)_ where Z₍𝐇_ℓ₎(_x_) = ∏ᵢ₌₀^(ℓ-1)(_x_ - ωₗⁱ).

1. **Prover:** Compute quotient _q(x)_ = _f(x)_/Z₍𝐇_ℓ₎(_x_) and send _Cᵩ_ = PCS(_q_)
2. **Verifier:** Send random challenge _r_ ∈ 𝔽_p
3. **Prover:** Open _f(r)_ and _q(r)_ with respective proofs
4. **Verifier:** Check _f(r)_ = _q(r)_ · Z₍𝐇_ℓ₎(_r_) and verify opening proofs

### Complexity Analysis

- **Prover time:** _O(D)_𝔾 + _O(D)_𝔽
- **Verifier time:** _O(1)_𝔾 + _O(1)_𝔽
- **Proof size:** _O(1)_

## Univariate SumCheck PIOP

A PIOP proving that polynomial evaluations sum to zero over subgroup 𝐇_ℓ.

### Relation

<img src="https://latex.codecogs.com/svg.latex?\mathcal{R}_{\text{UniSC}}%20=%20\left\{(C%20\in%20\mathbb{G};%20f(X)%20\in%20\mathbb{F}_p[X])%20:%20C%20=%20\text{PCS}(f)%20\land%20\sum_{i=0}^{\ell-1}%20f(\omega_\ell^i)%20=%200\right\}" />

### Protocol Description

**Key insight:** For polynomial _g(x)_ of degree < ℓ: ∑_(a ∈ 𝐇_ℓ) _g(a)_ = _g(0)_ · ℓ.

Therefore, ∑_(a ∈ 𝐇_ℓ) _f(a)_ = 0 ⟺ ∃ _h*_(x), _t(x)_ such that:

<img src="https://latex.codecogs.com/svg.latex?f(x)%20=%20h^*(x)%20\cdot%20Z_{\mathbb{H}_\ell}(x)%20+%20x%20\cdot%20t(x)" />

where deg(_t_) < ℓ - 1.

1. **Prover:** Compute _h*_(x) and _t(x)_, send _C_{h*}_ = PCS(_h*_) and _Cₜ_ = PCS(_t_)
2. **Verifier:** Send random challenge _r_ ∈ 𝔽_p
3. **Prover:** Open _f(r)_, _h*_(r), and _t(r)_ with respective proofs
4. **Verifier:** Check _f(r)_ = _h*_(r) · Z₍𝐇_ℓ₎(_r_) + _r_ · _t(r)_ and verify opening proofs

### Complexity Analysis

- **Prover time:** _O(D)_𝔾 + _O(D)_𝔽
- **Verifier time:** _O(1)_𝔾 + _O(1)_𝔽
- **Proof size:** _O(1)_

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

- **Library:** Used mcl library with BN_SNARK1 curve for cryptographic operations
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