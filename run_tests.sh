#!/bin/bash
cd build
export LD_LIBRARY_PATH=../mcl/lib:$LD_LIBRARY_PATH

echo "Running test suite..."
./test_suite

echo ""
echo "Running main program..."
./crypto_practice
