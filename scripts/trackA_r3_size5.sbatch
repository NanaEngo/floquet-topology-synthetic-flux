#!/bin/bash
#SBATCH --job-name=trackA_r3_size5
#SBATCH --partition=production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r3_size5_%j.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r3_size5_%j.err

# Audit run R3 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 1): fifth size,
# 5 cells / 10 modes (dim 3^10 = 59049), N_MC=300, disorders {0, 0.5},
# both loss profiles, both sectors, production seed 20260814.  The size-5
# contrast is scored OUT OF SAMPLE against data/r3_size5_prediction.json,
# which was frozen before this submission (L3 of the refinement plan).
# Terminal 3x3 reduced blocks are dumped so the R8 witness analysis can be
# extended to 5 cells without another dynamics run.
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
mkdir -p "$ROOT/logs" "$ROOT/data/chiral_production"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$PY" "$ROOT/scripts/run_r3_r8_records.py" \
  --n-cells 5 \
  --n-mc 300 \
  --disorder-levels 0.0 0.5 \
  --dump-blocks \
  --n-workers 16 \
  --seed 20260814 \
  --output "$ROOT/data/chiral_production/chiral_floquet_5cells_NMC300.json"
