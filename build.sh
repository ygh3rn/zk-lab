#!/bin/bash

# Cryptographic Programming Practice Build Script
# This script automates the setup and compilation process

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}Cryptographic Programming Practice Build Script${NC}"
echo "================================================="

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if we're in the right directory
if [[ ! -f "src/main.cpp" ]]; then
    print_error "Please run this script from the project root directory"
    print_error "Expected directory structure:"
    print_error "  crypto_practice/"
    print_error "    ├── src/"
    print_error "    ├── include/"
    print_error "    ├── mcl/"
    print_error "    └── build.sh"
    exit 1
fi

# Check for required dependencies
print_status "Checking dependencies..."

check_command() {
    if command -v $1 &> /dev/null; then
        print_status "$1 is installed"
        return 0
    else
        print_error "$1 is not installed"
        return 1
    fi
}

# Check for required tools
MISSING_DEPS=false
check_command "g++" || MISSING_DEPS=true
check_command "cmake" || MISSING_DEPS=true
check_command "make" || MISSING_DEPS=true
check_command "git" || MISSING_DEPS=true

if $MISSING_DEPS; then
    print_error "Please install missing dependencies:"
    echo "  Ubuntu/Debian: sudo apt install build-essential cmake git libgmp-dev"
    echo "  macOS: brew install cmake git gmp"
    exit 1
fi

# Setup MCL library if not already present
if [[ ! -d "mcl" ]]; then
    print_status "MCL library not found. Cloning from repository..."
    git clone https://github.com/herumi/mcl.git
    
    print_status "Building MCL library..."
    cd mcl
    
    # Build MCL with optimizations for BN_SNARK1
    make -j$(nproc) MCL_USE_LLVM=1 2>/dev/null || make -j$(nproc)
    
    # Run MCL tests to verify installation
    print_status "Testing MCL installation..."
    if ./bin/bn_test > /dev/null 2>&1; then
        print_status "MCL library built and tested successfully"
    else
        print_warning "MCL tests failed, but continuing with build"
    fi
    
    cd ..
else
    print_status "MCL library found"
fi

# Create build directory
print_status "Setting up build directory..."
mkdir -p build
cd build

# Configure with CMake
print_status "Configuring build with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build the project
print_status "Building project..."
make -j$(nproc)

# Check if build was successful
if [[ -f "crypto_practice" ]]; then
    print_status "Build completed successfully!"
    
    # Run basic test
    print_status "Running basic functionality test..."
    if ./crypto_practice > /dev/null 2>&1; then
        print_status "Basic functionality test passed"
    else
        print_warning "Basic functionality test failed - check implementation"
    fi
    
    # Build test suite if available
    if [[ -f "../test_suite.cpp" ]]; then
        print_status "Building test suite..."
        g++ -std=c++17 -O3 -DNDEBUG -DMCL_FP_BIT=384 -DMCL_FR_BIT=256 \
            -I../include -I../mcl/include \
            ../test_suite.cpp ../src/ntt.cpp ../src/polynomial.cpp ../src/kzg.cpp ../src/piop.cpp \
            -L../mcl/lib -lmcl -lgmp -lgmpxx -pthread \
            -o test_suite
        
        if [[ -f "test_suite" ]]; then
            print_status "Test suite built successfully"
            print_status "Running comprehensive tests..."
            LD_LIBRARY_PATH=../mcl/lib:$LD_LIBRARY_PATH ./test_suite
        fi
    fi
    
else
    print_error "Build failed!"
    exit 1
fi

cd ..

# Display usage information
echo ""
echo -e "${BLUE}Build completed successfully!${NC}"
echo "Usage:"
echo "  cd build"
echo "  ./crypto_practice          # Run main program"
echo "  ./test_suite               # Run comprehensive tests (if available)"
echo ""
echo "Project structure:"
echo "  crypto_practice/"
echo "    ├── src/                 # Source files"
echo "    ├── include/             # Header files"
echo "    ├── mcl/                 # MCL cryptography library"
echo "    ├── build/               # Build artifacts"
echo "    └── CMakeLists.txt       # Build configuration"
echo ""
echo "Implementation includes:"
echo "  ✓ Non-recursive NTT and inverse NTT"
echo "  ✓ Polynomial multiplication using NTT"
echo "  ✓ KZG commitment protocol with setup, commit, and verify"
echo "  ✓ Univariate ZeroTest PIOP"
echo "  ✓ Univariate SumCheck PIOP"
echo "  ✓ Performance benchmarks and comprehensive tests"
echo ""
echo -e "${GREEN}Ready for code review and demonstration!${NC}"