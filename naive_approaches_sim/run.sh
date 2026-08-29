#!/usr/bin/env bash
# Usage: bash run.sh [number of cores]

set -euo pipefail
cd "$(dirname "$0")"

NCORES=${1:-15}
NTASKS=400

# Avoid nested numerical-library parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

echo "=== Step 1: Generate shared data ==="
Rscript naive_setup.R

mkdir -p results

# Naive approaches (a)-(c)
echo "=== Step 2: Run naive approaches (a-c) — ${NCORES} cores ==="
for task_id in $(seq 1 $NTASKS); do
  Rscript naive_worker.R "$task_id" > "results/log_naive_${task_id}.txt" 2>&1 &
  while (( $(jobs -rp | wc -l) >= NCORES )); do
    sleep 1
  done
done
wait
echo "All $NTASKS naive tasks complete."

echo "=== Step 3: Report naive results ==="
Rscript naive_collect.R naive

# R-native S(ME)
echo "=== Step 4: Run signed KF (d), R-native S(ME) — ${NCORES} cores ==="
for task_id in $(seq 1 $NTASKS); do
  Rscript signed_worker.R "$task_id" > "results/log_signed_${task_id}.txt" 2>&1 &
  while (( $(jobs -rp | wc -l) >= NCORES )); do
    sleep 1
  done
done
wait
echo "All $NTASKS signed tasks complete."

# RG(signed, ME)
echo "=== Step 5: Run RG-KF (e) — ${NCORES} cores ==="
for task_id in $(seq 1 $NTASKS); do
  Rscript rg_worker.R "$task_id" > "results/log_rg_${task_id}.txt" 2>&1 &
  while (( $(jobs -rp | wc -l) >= NCORES )); do
    sleep 1
  done
done
wait
echo "All $NTASKS RG tasks complete."

echo "=== Step 6: Collect all results ==="
Rscript naive_collect.R all

echo "=== Step 7: Generate plot ==="
Rscript plot_naive_approaches.R

echo "=== Done ==="
