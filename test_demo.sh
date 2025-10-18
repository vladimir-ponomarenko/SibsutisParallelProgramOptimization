#!/bin/bash

echo "=== DNA Motif Finder Demo ==="
echo

echo "1. Тест с 1 MPI процессом:"
mpirun -n 1 ./build/DNAMotifFinder data/FoxA2_5000.fst data/FoxA2_major_30.mot | head -10
echo

echo "2. Тест с 2 MPI процессами:"
mpirun -n 2 ./build/DNAMotifFinder data/FoxA2_5000.fst data/FoxA2_major_30.mot | head -10
echo

echo "3. Тест с OpenMP (4 потока):"
OMP_NUM_THREADS=4 mpirun -n 2 ./build/DNAMotifFinder data/FoxA2_5000.fst data/FoxA2_major_30.mot | head -10
echo

echo "4. Тест с сохранением в файл:"
mpirun -n 2 ./build/DNAMotifFinder data/FoxA2_5000.fst data/FoxA2_major_30.mot demo_results.txt
echo "Результаты сохранены в demo_results.txt"
echo

echo "=== Демонстрация завершена ==="
