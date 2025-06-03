cmake_minimum_required(VERSION 3.12)
project(CryptographyPractice)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)

# Add MCL library
# Assuming MCL is installed in /usr/local or via package manager
find_path(MCL_INCLUDE_DIR mcl/bn.h 
    PATHS /usr/local/include /usr/include
    PATH_SUFFIXES include)

find_library(MCL_LIBRARY mcl
    PATHS /usr/local/lib /usr/lib
    PATH_SUFFIXES lib)

if(NOT MCL_INCLUDE_DIR OR NOT MCL_LIBRARY)
    message(FATAL_ERROR "MCL library not found. Please install MCL library.")
endif()

# Include directories
include_directories(${MCL_INCLUDE_DIR})

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Source files
set(SOURCES
    protocol.cpp
    kzg.cpp
)

# Create executable
add_executable(cryptography_practice ${SOURCES})
add_executable(complexity_benchmark complexity_benchmark.cpp kzg.cpp)
target_link_libraries(complexity_benchmark ${MCL_LIBRARY})

# Link libraries
target_link_libraries(cryptography_practice ${MCL_LIBRARY})

# Optional: Add test target
enable_testing()
add_test(NAME cryptography_practice_test COMMAND cryptography_practice)

# Installation
install(TARGETS cryptography_practice DESTINATION bin)
install(FILES kzg.h DESTINATION include)

# Print configuration info
message(STATUS "MCL Include Dir: ${MCL_INCLUDE_DIR}")
message(STATUS "MCL Library: ${MCL_LIBRARY}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")

# Custom target for running the program
add_custom_target(run
    COMMAND cryptography_practice
    DEPENDS cryptography_practice
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)