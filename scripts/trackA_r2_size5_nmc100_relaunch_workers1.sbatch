#!/bin/bash
#SBATCH --job-name=trackA_r2_size5_nmc100_w1
#SBATCH --partition=production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r2_size5_nmc100_w1_\%j.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r2_size5_nmc100_w1_\%j.err

# Relaunch of the stalled R2 cross-solver contrast check, SINGLE-WORKER variant:
# 5 cells / 10 modes, n_mc=100 (200 records = 100 samples x 2 sectors),
# ntraj=200, production seed 20260814, edge_low_loss, disorder 0.0,
# gamma_edge=1e-4, gamma_bulk=0.2, n_max=2, kerr=0.5, wigner_grid=31,
# n_time=9, n_periods=1.  SINGLE WORKER (--workers 1) to avoid the
# multiprocessing futex/GIL stall that killed the 16-worker attempts.
# Checkpointable via --fragment-dir.
# Output: data/solver_contrast_check_5cells_n200_nmc100.json
# (preserves existing solver_contrast_check_5cells_n200.json, n_mc=8/ntraj=200).
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
mkdir -p "$ROOT/logs" "$ROOT/data"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$PY" "$ROOT/scripts/run_solver_contrast_check.py"   --n-cells 5   --n-mc 100   --ntraj 200   --workers 1   --profile edge_low_loss   --disorder 0.0   --gamma-edge 1e-4   --gamma-bulk 0.2   --n-max 2   --kerr 0.5   --wigner-grid 31   --n-time 9   --n-periods 1   --output "$ROOT/data/solver_contrast_check_5cells_n200_nmc100.json"   --fragment-dir "$ROOT/_fragments/solver_contrast_5cells_nmc100"
