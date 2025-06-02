# ZK-Poly-Toolkit

A cryptographic library implementing polynomial commitment schemes and zero-knowledge Interactive Oracle Proofs. This project demonstrates core concepts in modern cryptographic protocols used in blockchain and privacy-preserving systems.

## Features

* NTT-based polynomial operations with O(n log n) complexity
* KZG polynomial commitments with constant-size proofs
* Zero-knowledge PIOPs including ZeroTest and SumCheck protocols
* Comprehensive test suite with performance benchmarks
* Automated build system

## Installation

### Prerequisites

**Ubuntu/Debian:**
```bash
sudo apt install build-essential cmake git libgmp-dev pkg-config
```

**macOS:**
```bash
brew install cmake git gmp
```

### Build and Run

```bash
git clone https://github.com/ygh3rn/zk-poly-toolkit.git
cd zk-poly-toolkit
./build.sh
./run_tests.sh
```

The build script automatically downloads and compiles the required MCL library for pairing-based cryptography.

## Implementation Details

### Number Theoretic Transform (NTT)

Implements the Cooley-Tukey algorithm for efficient polynomial multiplication in finite fields. The non-recursive approach achieves O(n log n) complexity for polynomial operations.

### KZG Polynomial Commitments

Based on the Kate-Zaverucha-Goldberg scheme, providing:

* **Setup**: Generates structured reference strings
* **Commit**: Creates constant-size polynomial commitments
* **CreateWitness**: Produces evaluation proofs
* **VerifyEval**: Verifies proofs using pairing operations
* **Batch operations**: Efficient verification of multiple evaluations

### Polynomial Interactive Oracle Proofs

* **ZeroTest PIOP**: Proves that a polynomial evaluates to zero on a multiplicative subgroup
* **SumCheck PIOP**: Proves that the sum of polynomial evaluations equals a claimed value

Both protocols achieve the theoretical optimal complexity bounds.

## Performance

Benchmarks on modern hardware show the implementation meets expected performance characteristics:

| Operation | Complexity | Time (n=1024) |
|-----------|------------|---------------|
| NTT Round-trip | O(n log n) | ~0.5ms |
| KZG Commitment | O(n) | ~18ms |
| KZG Verification | O(1) | ~0.6ms |

## Technical Stack

* C++17 with modern practices
* MCL library for elliptic curve operations
* CMake build system
* GMP for arbitrary precision arithmetic

## Applications

This implementation can serve as a foundation for:

* Zero-knowledge proof systems (SNARKs/STARKs)
* Blockchain scalability research
* Educational cryptography projects
* Privacy-preserving protocol development

## Important Security Notice

This implementation is designed for educational and research purposes. The KZG trusted setup generates parameters locally, which is not suitable for production use. Real-world applications require multi-party trusted setup ceremonies to ensure security.

## Testing

The test suite covers all major components and includes both correctness and performance tests. Run with:

```bash
./run_tests.sh
```

Expected output shows all tests passing with detailed performance metrics.

## References

* Thaler, J. (2022). Proofs, Arguments, and Zero-Knowledge. Available at: [https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf](https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf)
* Kate, A., Zaverucha, G. M., & Goldberg, I. (2010). Constant-Size Commitments to Polynomials and Their Applications. ASIACRYPT 2010.

## Contributing

This project was developed as part of cryptographic protocol research. The code is available for educational use and further development.

## License

See LICENSE file for details.
