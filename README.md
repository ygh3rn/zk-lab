# Cryptography Programming Practice

This repository contains implementations of cryptographic protocols including NTT, KZG commitments, and PIOPs using the MCL library.

## 📋 Requirements

- C++17 compatible compiler (g++ or clang++)
- MCL library (https://github.com/herumi/mcl)
- GMP library
- OpenSSL (for random number generation)

## 🛠️ Installation

### Installing MCL Library

```bash
# Clone MCL repository
git clone https://github.com/herumi/mcl
cd mcl

# Build and install
make -j4
sudo make install
```

### Building the Project

```bash
# Clone this repository
cd crypto-practice

# Build all components
make all

# Or build individual components
make test_ntt    # NTT tests only
make test_kzg    # KZG tests only
make test_piop   # PIOP tests only
```

## 🚀 Running Tests

### Run All Tests
```bash
make test
# or
./test_all
```

### Run Individual Test Suites
```bash
./test_ntt    # Test NTT implementation
./test_kzg    # Test KZG protocol
./test_piop   # Test ZeroTest and SumCheck PIOPs
```

## 📁 Project Structure

```
crypto-practice/
├── crypto_practice.hpp    # Main header file with declarations
├── ntt.cpp               # NTT implementation
├── polynomial_ops.cpp    # Polynomial operations
├── kzg.cpp              # KZG commitment scheme
├── zerotest.cpp         # ZeroTest PIOP
├── sumcheck.cpp         # SumCheck PIOP
├── test_ntt.cpp         # NTT test suite
├── test_kzg.cpp         # KZG test suite
├── test_piop.cpp        # PIOP test suite
├── test_all.cpp         # Main test runner
├── Makefile             # Build configuration
└── README.md            # This file
```

## 📊 Implemented Features

### 1. **NTT (Number Theoretic Transform)**
- Non-recursive Cooley-Tukey algorithm
- Inverse NTT
- Polynomial interpolation
- Time complexity: O(n log n)

### 2. **Polynomial Operations**
- Multiplication using NTT
- Evaluation using Horner's method
- Division with remainder
- Time complexity: O(n log n) for multiplication

### 3. **KZG Commitment Scheme**
- Setup with trusted parameters
- Polynomial commitment
- Witness creation for evaluations
- Pairing-based verification
- Prover time: O(D) group operations
- Verifier time: O(1) pairings

### 4. **ZeroTest PIOP**
- Proves polynomial evaluates to zero on a subgroup
- Uses quotient polynomial technique
- Constant proof size
- Prover: O(D) operations
- Verifier: O(1) operations

### 5. **SumCheck PIOP**
- Proves sum of evaluations equals zero
- Constant proof size
- Prover: O(D) operations
- Verifier: O(1) operations

## 🧪 Testing

Each component includes comprehensive tests:

- **Correctness tests**: Verify mathematical properties
- **Edge cases**: Test with zero polynomials, constant polynomials
- **Random tests**: Use randomly generated inputs
- **Performance benchmarks**: Measure time complexity

## ⚡ Performance

The implementation achieves the required time complexities:

| Operation | Time Complexity | Space Complexity |
|-----------|----------------|------------------|
| NTT | O(n log n) | O(n) |
| Polynomial Multiply | O(n log n) | O(n) |
| KZG Commit | O(D) | O(1) |
| KZG Verify | O(1) | O(1) |
| PIOP Prove | O(D) | O(D) |
| PIOP Verify | O(1) | O(1) |

## 🔧 Troubleshooting

### MCL Library Not Found
```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH
```

### Compilation Errors
Ensure you have C++17 support:
```bash
g++ --version  # Should be 7.0 or higher
```

## 📚 References

1. KZG Paper: [Constant-Size Commitments to Polynomials and Their Applications](https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf)
2. MCL Library: [Documentation](https://github.com/herumi/mcl/blob/master/api.md)
3. Number Theoretic Transform: [FFT over Finite Fields](https://en.wikipedia.org/wiki/Discrete_Fourier_transform_(general))

## 📝 Notes

- The implementation uses the BN_SNARK1 curve as specified
- Random challenges in PIOPs are simulated using CSPRNG
- In production, use Fiat-Shamir transform for non-interactivity
- The trusted setup in KZG should use a secure ceremony in practice

## ✅ Completion Checklist

- [x] Non-recursive NTT and inverse NTT
- [x] Polynomial interpolation example
- [x] Polynomial multiplication using NTT
- [x] KZG protocol implementation
- [x] Univariate ZeroTest PIOP
- [x] Univariate SumCheck PIOP
- [x] Comprehensive test suite
- [x] Performance benchmarks
- [x] Documentation