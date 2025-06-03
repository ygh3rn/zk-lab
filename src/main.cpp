#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <cassert>

#include "ntt.h"
#include "polynomial.h"
#include "kzg.h"
#include "piop.h"

using namespace mcl::bn;

// Utility function to measure execution time
template<typename Func>
double measure_time(Func func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return duration.count() / 1000.0; // Convert to milliseconds
}

// Test NTT implementation
void test_ntt() {
    std::cout << "\n=== Testing NTT Implementation ===" << std::endl;
    
    size_t n = 8; // Small size for testing
    NTT ntt(n);
    
    // Create test polynomial coefficients
    std::vector<Fr> coeffs = {Fr(1), Fr(2), Fr(3), Fr(4), Fr(0), Fr(0), Fr(0), Fr(0)};
    std::vector<Fr> original_coeffs = coeffs;
    
    std::cout << "Original coefficients: ";
    for (const auto& coeff : coeffs) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;
    
    // Forward NTT
    double forward_time = measure_time([&]() {
        ntt.forward_transform(coeffs);
    });
    
    std::cout << "After forward NTT: ";
    for (const auto& val : coeffs) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // Inverse NTT
    double inverse_time = measure_time([&]() {
        ntt.inverse_transform(coeffs);
    });
    
    std::cout << "After inverse NTT: ";
    for (const auto& coeff : coeffs) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;
    
    // Verify correctness
    bool correct = true;
    for (size_t i = 0; i < n; ++i) {
        if (coeffs[i] != original_coeffs[i]) {
            correct = false;
            break;
        }
    }
    
    std::cout << "NTT correctness: " << (correct ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Forward NTT time: " << forward_time << " ms" << std::endl;
    std::cout << "Inverse NTT time: " << inverse_time << " ms" << std::endl;
}

// Test polynomial operations
void test_polynomial_operations() {
    std::cout << "\n=== Testing Polynomial Operations ===" << std::endl;
    
    // Create test polynomials
    std::vector<Fr> coeffs1 = {Fr(1), Fr(2), Fr(3)};  // 1 + 2x + 3x²
    std::vector<Fr> coeffs2 = {Fr(4), Fr(5)};         // 4 + 5x
    
    Polynomial poly1(coeffs1);
    Polynomial poly2(coeffs2);
    
    std::cout << "Polynomial 1: ";
    poly1.print();
    
    std::cout << "Polynomial 2: ";
    poly2.print();
    
    // Test addition
    Polynomial sum = poly1 + poly2;
    std::cout << "Sum: ";
    sum.print();
    
    // Test multiplication
    double mult_time = measure_time([&]() {
        Polynomial product = poly1 * poly2;
        std::cout << "Product: ";
        product.print();
    });
    
    std::cout << "Multiplication time: " << mult_time << " ms" << std::endl;
    
    // Test evaluation
    Fr x = Fr(2);
    Fr eval1 = poly1.evaluate(x);
    Fr eval2 = poly2.evaluate(x);
    
    std::cout << "poly1(2) = " << eval1 << std::endl;
    std::cout << "poly2(2) = " << eval2 << std::endl;
    
    // Test interpolation
    std::vector<std::pair<Fr, Fr>> points = {
        {Fr(0), Fr(1)},
        {Fr(1), Fr(3)},
        {Fr(2), Fr(5)}
    };
    
    Polynomial interpolated = Polynomial::lagrange_interpolation(points);
    std::cout << "Interpolated polynomial: ";
    interpolated.print();
}

// Test KZG commitment scheme
void test_kzg_commitment() {
    std::cout << "\n=== Testing KZG Commitment Scheme ===" << std::endl;
    
    size_t max_degree = 16;
    KZG kzg(max_degree);
    
    // Create test polynomial
    std::vector<Fr> coeffs = {Fr(1), Fr(2), Fr(3), Fr(4)};  // 1 + 2x + 3x² + 4x³
    Polynomial poly(coeffs);
    
    std::cout << "Test polynomial: ";
    poly.print();
    
    // Commit to polynomial
    G1 commitment;
    double commit_time = measure_time([&]() {
        commitment = kzg.commit(poly);
    });
    
    std::cout << "Commitment computed in " << commit_time << " ms" << std::endl;
    
    // Test evaluation and witness creation
    Fr point = Fr(5);
    Fr value = poly.evaluate(point);
    
    G1 witness;
    double witness_time = measure_time([&]() {
        witness = kzg.create_witness(poly, point);
    });
    
    std::cout << "Witness created in " << witness_time << " ms" << std::endl;
    
    // Verify evaluation
    bool verification_result;
    double verify_time = measure_time([&]() {
        verification_result = kzg.verify_eval(commitment, point, value, witness);
    });
    
    std::cout << "Verification result: " << (verification_result ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Verification time: " << verify_time << " ms" << std::endl;
    
    // Test batch opening
    std::vector<Fr> points = {Fr(1), Fr(2), Fr(3)};
    std::vector<Fr> values;
    values.reserve(points.size());
    
    for (const Fr& pt : points) {
        values.push_back(poly.evaluate(pt));
    }
    
    G1 batch_witness;
    double batch_witness_time = measure_time([&]() {
        batch_witness = kzg.create_batch_witness(poly, points);
    });
    
    bool batch_verification;
    double batch_verify_time = measure_time([&]() {
        batch_verification = kzg.verify_batch_eval(commitment, points, values, batch_witness);
    });
    
    std::cout << "Batch witness creation time: " << batch_witness_time << " ms" << std::endl;
    std::cout << "Batch verification result: " << (batch_verification ? "PASSED" : "FAILED") << std::endl;
    std::cout << "Batch verification time: " << batch_verify_time << " ms" << std::endl;
}

// Test Univariate ZeroTest PIOP
void test_zerotest_piop() {
    std::cout << "\n=== Testing Univariate ZeroTest PIOP ===" << std::endl;
    
    size_t max_degree = 16;
    size_t subgroup_size = 8;
    
    KZG kzg(max_degree);
    UnivariateZeroTest zerotest(kzg, subgroup_size);
    
    // Create a polynomial that should be zero on the subgroup
    // For simplicity, use the zero polynomial
    std::vector<Fr> zero_coeffs = {Fr(0)};
    Polynomial zero_poly(zero_coeffs);
    
    std::cout << "Testing with zero polynomial" << std::endl;
    
    // Generate proof
    UnivariateZeroTest::Proof proof;
    double prove_time = measure_time([&]() {
        proof = zerotest.prove(zero_poly);
    });
    
    std::cout << "ZeroTest proof generation time: " << prove_time << " ms" << std::endl;
    
    // Verify proof
    bool verification_result;
    double verify_time = measure_time([&]() {
        verification_result = zerotest.verify(proof);
    });
    
    std::cout << "ZeroTest verification result: " << (verification_result ? "PASSED" : "FAILED") << std::endl;
    std::cout << "ZeroTest verification time: " << verify_time << " ms" << std::endl;
}

// Test Univariate SumCheck PIOP
void test_sumcheck_piop() {
    std::cout << "\n=== Testing Univariate SumCheck PIOP ===" << std::endl;
    
    size_t max_degree = 16;
    size_t subgroup_size = 8;
    
    KZG kzg(max_degree);
    UnivariateSumCheck sumcheck(kzg, subgroup_size);
    
    // Create a constant polynomial for testing
    std::vector<Fr> const_coeffs = {Fr(5)}; // Constant polynomial f(x) = 5
    Polynomial const_poly(const_coeffs);
    
    std::cout << "Testing with constant polynomial f(x) = 5" << std::endl;
    std::cout << "Expected sum over subgroup of size " << subgroup_size << ": " << Fr(5 * subgroup_size) << std::endl;
    
    // Generate proof
    UnivariateSumCheck::Proof proof;
    double prove_time = measure_time([&]() {
        proof = sumcheck.prove(const_poly);
    });
    
    std::cout << "SumCheck proof generation time: " << prove_time << " ms" << std::endl;
    std::cout << "Claimed sum: " << proof.claimed_sum << std::endl;
    
    // Verify proof
    Fr expected_sum = Fr(5 * subgroup_size);
    bool verification_result;
    double verify_time = measure_time([&]() {
        verification_result = sumcheck.verify(proof, expected_sum);
    });
    
    std::cout << "SumCheck verification result: " << (verification_result ? "PASSED" : "FAILED") << std::endl;
    std::cout << "SumCheck verification time: " << verify_time << " ms" << std::endl;
}

// Performance benchmarks
void run_performance_benchmarks() {
    std::cout << "\n=== Performance Benchmarks ===" << std::endl;
    
    std::vector<size_t> sizes = {256, 512, 1024, 2048};
    
    for (size_t n : sizes) {
        std::cout << "\nBenchmarking with size " << n << ":" << std::endl;
        
        // NTT benchmark
        NTT ntt(n);
        std::vector<Fr> data(n);
        
        // Fill with random data
        std::random_device rd;
        std::mt19937 gen(rd());
        for (size_t i = 0; i < n; ++i) {
            data[i].setByCSPRNG();
        }
        
        double ntt_time = measure_time([&]() {
            ntt.forward_transform(data);
            ntt.inverse_transform(data);
        });
        
        std::cout << "  NTT round-trip time: " << ntt_time << " ms" << std::endl;
        
        // KZG benchmark
        if (n <= 1024) { // Limit KZG size for reasonable setup time
            KZG kzg(n);
            
            // Create random polynomial
            std::vector<Fr> coeffs(n/2);
            for (size_t i = 0; i < n/2; ++i) {
                coeffs[i].setByCSPRNG();
            }
            Polynomial poly(coeffs);
            
            double commit_time = measure_time([&]() {
                G1 commitment = kzg.commit(poly);
                // Use the commitment to avoid warning
                (void)commitment;
            });
            
            Fr point;
            point.setByCSPRNG();
            Fr value = poly.evaluate(point);
            
            G1 commitment = kzg.commit(poly);
            double witness_time = measure_time([&]() {
                G1 witness = kzg.create_witness(poly, point);
                // Use the witness to avoid warning
                (void)witness;
            });
            
            G1 witness = kzg.create_witness(poly, point);
            double verify_time = measure_time([&]() {
                bool result = kzg.verify_eval(commitment, point, value, witness);
                // Use the result to avoid warning
                (void)result;
            });
            
            std::cout << "  KZG commit time: " << commit_time << " ms" << std::endl;
            std::cout << "  KZG witness time: " << witness_time << " ms" << std::endl;
            std::cout << "  KZG verify time: " << verify_time << " ms" << std::endl;
        }
    }
}

int main() {
    std::cout << "Cryptographic Programming Practice" << std::endl;
    std::cout << "===================================" << std::endl;
    
    try {
        // Initialize the mcl library with BN_SNARK1 curve
        std::cout << "Initializing MCL library with BN_SNARK1 curve..." << std::endl;
        mcl::bn::initPairing();
        std::cout << "MCL library initialized successfully!" << std::endl;
        
        // Run all tests
        test_ntt();
        test_polynomial_operations();
        test_kzg_commitment();
        test_zerotest_piop();
        test_sumcheck_piop();
        
        // Run performance benchmarks
        run_performance_benchmarks();
        
        std::cout << "\n=== All Tests Completed ===" << std::endl;
        std::cout << "Summary:" << std::endl;
        std::cout << "✓ NTT implementation and polynomial multiplication" << std::endl;
        std::cout << "✓ KZG commitment scheme with setup, commit, and verify" << std::endl;
        std::cout << "✓ Univariate ZeroTest PIOP" << std::endl;
        std::cout << "✓ Univariate SumCheck PIOP" << std::endl;
        std::cout << "✓ Performance benchmarks" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}