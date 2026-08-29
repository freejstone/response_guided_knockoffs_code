#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <n> <k> <value>" >&2
  exit 1
fi

n="$1"
k="$2"
value="$3"

# Avoid oversubscribing BLAS/OpenMP when running many jobs
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

# Optional per-job logging (uncomment if desired)
# mkdir -p logs
# log="logs/val_${value}_n_${n}_A_${A}.out"
# exec >"$log" 2>&1

Rscript run_sim.R "$value" "$n" "$k"
