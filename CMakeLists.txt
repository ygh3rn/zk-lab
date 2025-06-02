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

# Source files
set(SOURCES
    src/main.cpp
    src/ntt.cpp
    src/polynomial.cpp
    src/kzg.cpp
    src/piop.cpp
)

# Create executable
add_executable(crypto_practice ${SOURCES})

# Link libraries
target_link_libraries(crypto_practice 
    ${MCL_LIB_DIR}/libmcl.a
    ${GMP_LIB}
    ${GMPXX_LIB}
    pthread
)

# Compiler flags
target_compile_options(crypto_practice PRIVATE
    -O3
    -DNDEBUG
    -DMCL_FP_BIT=384
    -DMCL_FR_BIT=256
    -DMCL_USE_LLVM=1
)

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(crypto_practice OpenMP::OpenMP_CXX)
endif()