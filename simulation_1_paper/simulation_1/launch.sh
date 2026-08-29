#!/usr/bin/env bash
set -euo pipefail

# Ensure we're in bash (not sh/zsh)
if [ -z "${BASH_VERSION:-}" ]; then
  echo "Please run with bash: bash $0 ..." >&2
  exit 1
fi

# Defaults (override via flags)
ns=(1005 1250 1500 1750 2000)
ks=(50 100 150)

# values default: 0..99 (Bash-3.2 friendly; no mapfile/readarray)
values=()
for i in $(seq 0 99); do values+=("$i"); done

jobs=20   # up to 20 concurrent jobs, change with -j

usage() {
  cat >&2 <<'USAGE'
Usage: launch.sh [--ns "0.3,0.7"] [--ks "3.5,7"]
                 [--values "0:99" | "0,5,10"] [-j 20]
Examples:
  ./launch.sh -j 20
  ./launch.sh --ns "0.3,0.7" --ks "3.5,7" --values "0:19" -j 20
  ./launch.sh --values "0,5,10,42" -j 8
USAGE
  exit 1
}

# Parse flags (Bash-3.2 safe)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ns)
      IFS=',' read -r -a ns <<< "$2"
      shift 2 ;;
    --ks|--as)
      IFS=',' read -r -a ks <<< "$2"
      shift 2 ;;
    --values)
      if [[ "$2" == *:* ]]; then
        IFS=':' read -r start end <<< "$2"
        values=()
        # Bash-3.2 doesn't support {start..end} with vars; use a loop
        i="$start"
        while [ "$i" -le "$end" ]; do
          values+=("$i")
          i=$((i+1))
        done
      else
        IFS=',' read -r -a values <<< "$2"
      fi
      shift 2 ;;
    -j|--jobs)
      jobs="$2"
      shift 2 ;;
    -h|--help)
      usage ;;
    *)
      echo "Unknown option: $1" >&2
      usage ;;
  esac
done

# Keep R's internal threading from oversubscribing when -j is large
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

triples=()
for n in "${ns[@]}"; do
  for k in "${ks[@]}"; do
    for v in "${values[@]}"; do
      triples+=("$n $k $v")
    done
  done
done
printf "%s\n" "${triples[@]}" | xargs -n3 -P "$jobs" ./run_one.sh
