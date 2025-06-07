cmake_minimum_required(VERSION 3.12)
project(CryptographyPractice)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_path(MCL_INCLUDE_DIR mcl/bn.hpp 
    PATHS /usr/local/include /usr/include
    PATH_SUFFIXES include)

find_library(MCL_LIBRARY mcl
    PATHS /usr/local/lib /usr/lib
    PATH_SUFFIXES lib)

if(NOT MCL_INCLUDE_DIR OR NOT MCL_LIBRARY)
    message(FATAL_ERROR "MCL library not found. Please install MCL first.")
endif()

include_directories(${MCL_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -march=native")

set(SOURCES
    src/ntt.cpp
    src/polynomial.cpp
    src/kzg.cpp
    src/zerotest.cpp
    src/sumcheck.cpp
)

add_executable(test_suite tests/test_suite.cpp ${SOURCES})
target_link_libraries(test_suite ${MCL_LIBRARY})

set_target_properties(test_suite PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

enable_testing()
add_test(NAME cryptography_test COMMAND test_suite)

add_custom_target(run
    COMMAND test_suite
    DEPENDS test_suite
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running cryptography test suite"
)

message(STATUS "MCL include directory: ${MCL_INCLUDE_DIR}")
message(STATUS "MCL library: ${MCL_LIBRARY}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CXX flags: ${CMAKE_CXX_FLAGS}")