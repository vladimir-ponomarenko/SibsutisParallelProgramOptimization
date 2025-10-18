#!/bin/bash

# This script builds the project and runs tests

set -e

echo "=== DNA Motif Finder Build Script ==="

if [ ! -f "CMakeLists.txt" ]; then
    echo "Error: CMakeLists.txt not found. Please run this script from the project root."
    exit 1
fi

echo "Creating build directory..."
mkdir -p build
cd build

echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

echo "Building project..."
make -j$(nproc)

echo "Build completed successfully!"

echo "Running tests..."
make test

echo "All tests passed!"

echo "Building examples..."
if [ -d "../examples" ]; then
    echo "Examples directory found."
else
    echo "Examples directory not found."
fi

echo "=== Build script completed successfully ==="
echo ""
echo "To run the application:"
echo "  mpirun -n 4 ./DNAMotifFinder data/FoxA2_5000.fst data/FoxA2_major_30.mot"
echo ""
echo "To run with custom settings:"
echo "  mpirun -n 4 ./DNAMotifFinder --threads 8 --verbose data/FoxA2_5000.fst data/FoxA2_major_30.mot results.txt"
