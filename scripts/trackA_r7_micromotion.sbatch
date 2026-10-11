#!/bin/bash
#SBATCH --job-name=trackA_r7_micro
#SBATCH --partition=production
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r7_micro_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r7_micro_%A_%a.err

# Audit run R7 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 1):
# sub-period (micromotion) scan: the observables are recorded at
# t = 0, T/8, ..., T instead of only at t = T, for both loss profiles and
# both disorder levels, at two and three cells.  Answers whether the sign and
# magnitude of the local contrast depend on when in the period the observable
# is read out (Xu2022 micromotion argument).
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

CELLS=$(( SLURM_ARRAY_TASK_ID + 2 ))

exec "$PY" "$ROOT/scripts/run_revision_scans.py" micromotion \
  --n-cells "$CELLS" \
  --n-mc 500 \
  --ntraj 10 \
  --disorder-levels 0.0 0.35 \
  --seed 20260817 \
  --output-dir "$ROOT/data/revision_scans"
