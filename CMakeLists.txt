cmake_minimum_required(VERSION 3.12)
project(CryptographyPractice)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)

find_path(MCL_INCLUDE_DIR mcl/bn.h 
    PATHS /usr/local/include /usr/include
    PATH_SUFFIXES include)

find_library(MCL_LIBRARY mcl
    PATHS /usr/local/lib /usr/lib
    PATH_SUFFIXES lib)

if(NOT MCL_INCLUDE_DIR OR NOT MCL_LIBRARY)
    message(FATAL_ERROR "MCL library not found")
endif()

include_directories(${MCL_INCLUDE_DIR} include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -march=native")

set(SOURCES
    src/ntt.cpp
    src/polynomial.cpp
    src/kzg.cpp
    src/piop.cpp
)

add_executable(test_suite tests/test_suite.cpp ${SOURCES})
target_link_libraries(test_suite ${MCL_LIBRARY})

enable_testing()
add_test(NAME cryptography_test COMMAND test_suite)

add_custom_target(run
    COMMAND test_suite
    DEPENDS test_suite
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR})