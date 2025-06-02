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

# Common source files
set(COMMON_SOURCES
    src/ntt.cpp
    src/polynomial.cpp
    src/kzg.cpp
    src/piop.cpp
)

# Main executable
add_executable(crypto_practice 
    src/main.cpp
    ${COMMON_SOURCES}
)

# Comprehensive test suite executable
add_executable(comprehensive_test_suite 
    test_comprehensive.cpp
    ${COMMON_SOURCES}
)

# Original test suite executable (if exists)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test_suite.cpp")
    add_executable(test_suite 
        test_suite.cpp
        ${COMMON_SOURCES}
    )
    
    # Link libraries for test_suite
    target_link_libraries(test_suite 
        ${MCL_LIB_DIR}/libmcl.a
        ${GMP_LIB}
        ${GMPXX_LIB}
        pthread
    )
    
    # Compiler flags for test_suite
    target_compile_options(test_suite PRIVATE
        -O3
        -DNDEBUG
        -DMCL_FP_BIT=384
        -DMCL_FR_BIT=256
        -DMCL_USE_LLVM=1
        -Wall
        -Wextra
    )
endif()

# Link libraries for main executable
target_link_libraries(crypto_practice 
    ${MCL_LIB_DIR}/libmcl.a
    ${GMP_LIB}
    ${GMPXX_LIB}
    pthread
)

# Link libraries for comprehensive test suite
target_link_libraries(comprehensive_test_suite 
    ${MCL_LIB_DIR}/libmcl.a
    ${GMP_LIB}
    ${GMPXX_LIB}
    pthread
)

# Compiler flags for main executable
target_compile_options(crypto_practice PRIVATE
    -O3
    -DNDEBUG
    -DMCL_FP_BIT=384
    -DMCL_FR_BIT=256
    -DMCL_USE_LLVM=1
    -Wall
    -Wextra
)

# Compiler flags for comprehensive test suite
target_compile_options(comprehensive_test_suite PRIVATE
    -O3
    -DNDEBUG
    -DMCL_FP_BIT=384
    -DMCL_FR_BIT=256
    -DMCL_USE_LLVM=1
    -Wall
    -Wextra
)

# Enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(crypto_practice OpenMP::OpenMP_CXX)
    target_link_libraries(comprehensive_test_suite OpenMP::OpenMP_CXX)
    if(TARGET test_suite)
        target_link_libraries(test_suite OpenMP::OpenMP_CXX)
    endif()
endif()

# Add custom target for running comprehensive tests
add_custom_target(run_comprehensive_tests
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/comprehensive_test_suite
    DEPENDS comprehensive_test_suite
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running comprehensive test suite with fixed implementations"
)

# Add custom target for running all tests
add_custom_target(run_all_tests
    COMMAND echo "Running comprehensive test suite..."
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/comprehensive_test_suite
    COMMAND echo "Running main program..."
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/crypto_practice
    DEPENDS comprehensive_test_suite crypto_practice
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running all tests and main program"
)