#!/bin/bash

# Simple test script for DNA Motif Finder

echo "=== DNA Motif Finder Simple Test ==="

# Test 1: Single process
echo "Test 1: Single process"
cd build
mpirun -n 1 ./DNAMotifFinder ../data/FoxA2_5000.fst ../data/FoxA2_major_30.mot > test1_output.txt
echo "Results saved to test1_output.txt"

# Test 2: Multiple processes
echo "Test 2: Multiple processes"
mpirun -n 2 ./DNAMotifFinder ../data/FoxA2_5000.fst ../data/FoxA2_major_30.mot test2_output.txt
echo "Results saved to test2_output.txt"

# Test 3: With OpenMP threads
echo "Test 3: With OpenMP threads"
mpirun -n 1 ./DNAMotifFinder --threads 4 --verbose ../data/FoxA2_5000.fst ../data/FoxA2_major_30.mot test3_output.txt
echo "Results saved to test3_output.txt"

echo "=== Test completed ==="

