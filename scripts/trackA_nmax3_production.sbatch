#!/bin/bash
#SBATCH --job-name=trackA_nmax3_prod
#SBATCH --partition=production
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_nmax3_prod_%j.log
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_nmax3_prod_%j.err

# Cutoff-convergence production: cells {2,3,4} at N_max=3, N_MC=500, same seeds
# as the N_max=2 production (nested matched design). Row-level checkpointing via
# run_8mode_nmc500_parallel.py (.partial.jsonl resume on requeue); finalized
# sizes are skipped. Exact Q<=2 subspace dynamics (validated vs full space
# 2026-10-02) — see the module docstring.
set -euo pipefail
# One BLAS thread per worker (avoids 24x24 oversubscription on the shared node)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
mkdir -p "$ROOT/logs" "$ROOT/data/chiral_production_nmax3"

for CELLS in 2 3 4; do
  echo "=== cells=$CELLS N_max=3 start $(date -Is)"
  "$PY" "$ROOT/scripts/run_8mode_nmc500_parallel.py" \
    --n-cells "$CELLS" \
    --n-max 3 \
    --n-mc 500 \
    --ntraj 10 \
    --n-workers 24 \
    --wigner-grid 31 \
    --kerr 0.5 \
    --gamma-edge 1e-4 \
    --gamma-bulk 0.2 \
    --t-stop 6.283185307179586 \
    --n-time 9 \
    --seed 20260814 \
    --output "$ROOT/data/chiral_production_nmax3/chiral_floquet_${CELLS}cells_NMC500_nmax3.json"
done
echo "NMAX3_PRODUCTION_DONE rc=$?"
