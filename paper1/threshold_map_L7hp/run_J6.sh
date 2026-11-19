#!/bin/zsh
# J = 6 thresholds with double-double assembly (2026-11-19, laptop). Scales reused from run 5 (paper1/threshold_map/data).
cd "$(dirname "$0")/../.."
export HP_NPROC=8 THREEBODY_BLAS_THREADS=8 HP_DATA_DIR=paper1/threshold_map_L7hp/data_J6
for T in J6v0 J6v1; do
  python3 paper1/threshold_map_L7hp/run_map_hp.py --target $T --stage 2 3 4 5
  echo "$T EXIT $?"
done
