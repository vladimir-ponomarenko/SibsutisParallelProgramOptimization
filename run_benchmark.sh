#!/bin/bash

# Disable warnings on Wayland
unset DISPLAY

EXE="./build/DNAMotifFinder"
DATA_SEQ="data/FoxA2_5000.fst"
DATA_MOTIF_ORIG="data/FoxA2_major_30.mot"
DATA_MOTIF_HEAVY="data/benchmark_extended.mot"

CSV_ORIG="results_original.csv"
CSV_HEAVY="results_extended.csv"
IMG_ORIG="graph_original.png"
IMG_HEAVY="graph_extended.png"

if [ ! -f "$EXE" ]; then echo "Ошибка: $EXE не найден."; exit 1; fi

run_suite() {
    local MOTIF_FILE=$1
    local OUT_FILE=$2
    local DESC=$3

    echo "----------------------------------------------------------------"
    echo "Тест: $DESC"
    echo "----------------------------------------------------------------"
    echo "Type,Cores,Time_Seconds" > "$OUT_FILE"

    # 1. MPI Scaling
    for np in 1 2 4 8 16; do
        output=$(OMP_NUM_THREADS=1 mpirun --oversubscribe --bind-to none -n $np $EXE $DATA_SEQ $MOTIF_FILE --verbose 2>/dev/null)
        time_val=$(echo "$output" | grep "total_processing_time" | awk '{print $2}' | head -n 1)
        if [ -n "$time_val" ]; then
            printf "MPI       | %-3d procs | %s s\n" "$np" "$time_val"
            echo "MPI,$np,$time_val" >> "$OUT_FILE"
        fi
    done

    # 2. OpenMP Scaling
    for thr in 1 2 4 8 16; do
        output=$(OMP_NUM_THREADS=$thr mpirun --oversubscribe --bind-to none -n 1 $EXE $DATA_SEQ $MOTIF_FILE --verbose --threads $thr 2>/dev/null)
        time_val=$(echo "$output" | grep "total_processing_time" | awk '{print $2}' | head -n 1)
        if [ -n "$time_val" ]; then
            printf "OpenMP    | %-3d thr   | %s s\n" "$thr" "$time_val"
            echo "OpenMP,$thr,$time_val" >> "$OUT_FILE"
        fi
    done

    # 3. Hybrid (2 Proc x N Threads)
    for thr in 1 2 4 8 16; do
        output=$(OMP_NUM_THREADS=$thr mpirun --oversubscribe --bind-to none -n 2 $EXE $DATA_SEQ $MOTIF_FILE --verbose --threads $thr 2>/dev/null)
        time_val=$(echo "$output" | grep "total_processing_time" | awk '{print $2}' | head -n 1)
        if [ -n "$time_val" ]; then
            printf "Hybrid 2x | %-3d thr   | %s s\n" "$thr" "$time_val"
            echo "Hybrid2,$thr,$time_val" >> "$OUT_FILE"
        fi
    done

    # 4. Hybrid (4 Proc x N Threads)
    for thr in 1 2 4 8 16; do
        output=$(OMP_NUM_THREADS=$thr mpirun --oversubscribe --bind-to none -n 4 $EXE $DATA_SEQ $MOTIF_FILE --verbose --threads $thr 2>/dev/null)
        time_val=$(echo "$output" | grep "total_processing_time" | awk '{print $2}' | head -n 1)
        if [ -n "$time_val" ]; then
            printf "Hybrid 4x | %-3d thr   | %s s\n" "$thr" "$time_val"
            echo "Hybrid4,$thr,$time_val" >> "$OUT_FILE"
        fi
    done
    echo ""
}

run_suite "$DATA_MOTIF_ORIG" "$CSV_ORIG" "Оригинальные данные"

echo "Генерация x150 данных..."
grep -v "^#" "$DATA_MOTIF_ORIG" > .temp_motifs
rm -f "$DATA_MOTIF_HEAVY"
for i in {1..150}; do cat .temp_motifs >> "$DATA_MOTIF_HEAVY"; done
rm .temp_motifs

run_suite "$DATA_MOTIF_HEAVY" "$CSV_HEAVY" "Расширенные данные"

if command -v python3 &>/dev/null; then
    python3 plot_graphs.py "$CSV_ORIG" "$IMG_ORIG" "Оригинальные данные"
    python3 plot_graphs.py "$CSV_HEAVY" "$IMG_HEAVY" "Расширенные данные (x150)"
fi

rm -f "$DATA_MOTIF_HEAVY"
echo "Готово."