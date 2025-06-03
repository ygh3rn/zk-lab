cmake_minimum_required(VERSION 3.12)
project(CryptographicPractice)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)
find_library(GMP_LIB gmp REQUIRED)
find_library(GMPXX_LIB gmpxx REQUIRED)

# MCL library path (adjust as needed)
set(MCL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/mcl")
set(MCL_INCLUDE_DIR "${MCL_ROOT}/include")
set(MCL_LIB_DIR "${MCL_ROOT}/lib")

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${MCL_INCLUDE_DIR})

# Common source files (excluding main.cpp)
set(COMMON_SOURCES
    src/ntt.cpp
    src/polynomial.cpp
    src/kzg.cpp
    src/piop.cpp
)

# Main executable
add_executable(crypto_practice ${COMMON_SOURCES} src/main.cpp)

# Test suite executable
add_executable(test_suite ${COMMON_SOURCES} test_suite.cpp)

# Link libraries for both executables
foreach(target crypto_practice test_suite)
    target_link_libraries(${target} 
        ${MCL_LIB_DIR}/libmcl.a
        ${GMP_LIB}
        ${GMPXX_LIB}
        pthread
    )
    
    # Compiler flags
    target_compile_options(${target} PRIVATE
        -O3
        -DNDEBUG
        -DMCL_FP_BIT=384
        -DMCL_FR_BIT=256
        -DMCL_USE_LLVM=1
        -Wall
        -Wextra
        -Wno-unused-parameter
    )
endforeach()

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(crypto_practice OpenMP::OpenMP_CXX)
    target_link_libraries(test_suite OpenMP::OpenMP_CXX)
endif()

# Create a custom target for comprehensive tests
add_custom_target(comprehensive_test_suite ALL
    DEPENDS test_suite
    COMMENT "Building comprehensive test suite"
)