# Cryptography Programming Practice

Complete implementation of cryptographic protocols including NTT, KZG commitments, and polynomial IOPs.

## Project Structure

```
cryptography_practice/
├── kzg.h                 # KZG protocol header
├── kzg.cpp              # KZG implementation
├── protocol.cpp         # Main implementation with NTT and PIOPs
├── CMakeLists.txt       # Build configuration
└── README.md            # This file
```

## Features Implemented

### 1. Number Theoretic Transform (NTT)
- Non-recursive NTT and inverse NTT
- Polynomial interpolation using inverse NTT
- Time complexity: O(n log n)

### 2. Polynomial Operations
- Polynomial multiplication using NTT
- Efficient polynomial division
- Polynomial evaluation using Horner's method

### 3. KZG Polynomial Commitment Scheme
- **Setup**: Generates structured reference string
- **Commit**: Creates polynomial commitments
- **CreateWitness**: Generates evaluation proofs
- **VerifyEval**: Verifies evaluation proofs using pairings

### 4. Univariate ZeroTest PIOP
Proves that a polynomial evaluates to zero on all points of a subgroup.
- **Prover time**: O(D)𝔾 + O(D)𝔽
- **Verifier time**: O(1)𝔾 + O(1)𝔽
- **Proof size**: O(1)

### 5. Univariate SumCheck PIOP
Proves that the sum of polynomial evaluations on a subgroup equals zero.
- **Prover time**: O(D)𝔾 + O(D)𝔽
- **Verifier time**: O(1)𝔾 + O(1)𝔽
- **Proof size**: O(1)

## Dependencies

### MCL Library
Install the MCL cryptographic library:

```bash
# Clone and build MCL
git clone https://github.com/herumi/mcl.git
cd mcl
make -j$(nproc)
sudo make install
```

### System Requirements
- C++17 compiler (GCC 7+ or Clang 6+)
- CMake 3.12+
- Linux/macOS (Windows with WSL)

## Compilation

```bash
# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the project
make -j$(nproc)

# Run the implementation
./cryptography_practice
```

## Usage

The program runs comprehensive tests for all implemented components:

```bash
./cryptography_practice
```

### Expected Output
The program will display:
1. NTT/INTT correctness verification and timing
2. Polynomial multiplication benchmarks
3. KZG protocol performance metrics
4. ZeroTest PIOP proof generation and verification
5. SumCheck PIOP proof generation and verification
6. Complete performance summary

## Performance Characteristics

### Time Complexity
- **NTT/INTT**: O(n log n) for size n
- **Polynomial Multiplication**: O(n log n)
- **KZG Setup**: O(n) group operations
- **KZG Commit**: O(n) group operations
- **KZG Prove**: O(n) field operations
- **KZG Verify**: O(1) pairing operations

### Proof Sizes
- KZG proof: ~80 bytes (1 G1 element + 1 Fr element)
- ZeroTest proof: ~80 bytes
- SumCheck proof: ~80 bytes

## Technical Details

### Curve Configuration
- Uses BN_SNARK1 curve as specified
- 254-bit prime field
- 128-bit security level
- Optimal ate pairing

### Implementation Notes
- All polynomial operations use coefficient representation
- NTT requires power-of-2 sizes for optimal performance
- Primitive roots of unity are computed for the BN_SNARK1 field
- Pairing verification uses optimized MCL library functions

## Testing

The implementation includes comprehensive tests:

1. **Unit Tests**: Individual component verification
2. **Integration Tests**: Full protocol execution
3. **Performance Tests**: Timing measurements
4. **Correctness Tests**: Mathematical property verification

### Random Testing
All tests use cryptographically secure random polynomials to ensure robustness.

## Security Considerations

- **Trusted Setup**: KZG requires trusted setup with toxic waste disposal
- **Field Choice**: BN_SNARK1 provides adequate security for research
- **Implementation**: Basic implementation for educational purposes

## References

1. Kate, A., Zaverucha, G. M., & Goldberg, I. (2010). Constant-size commitments to polynomials and their applications.
2. MCL Library: https://github.com/herumi/mcl
3. Number Theoretic Transform implementations
4. Polynomial IOPs and commitment schemes

## License

Educational implementation for cryptography practice assignment.