#!/bin/bash
#SBATCH --job-name=trackA_r9_longtime
#SBATCH --partition=production
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r9_longtime_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r9_longtime_%A_%a.err

# Audit run R9 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 1): stroboscopic
# negativity contrast curves t = 0..64T (production propagator, seeds
# 20260814) at 2 and 3 cells over all four disorder levels, N_MC=200, both
# profiles and sectors.  Determines whether the production sign change at
# t = T survives on the quasi-stationary plateau and on the way to the
# stationary state (referee question behind R9).
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
mkdir -p "$ROOT/logs" "$ROOT/data/r9_longtime"

case "${SLURM_ARRAY_TASK_ID}" in
  0) CELLS=2 ;;
  1) CELLS=3 ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$PY" "$ROOT/scripts/run_r9_longtime.py" \
  --n-cells "$CELLS" \
  --n-mc 200 \
  --horizon-periods 64 \
  --n-workers 16 \
  --seed 20260814 \
  --output "$ROOT/data/r9_longtime/r9_longtime_${CELLS}cells_NMC200.json"
