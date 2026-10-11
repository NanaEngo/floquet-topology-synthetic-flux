#!/bin/bash
#SBATCH --job-name=trackA_r2_xsolver
#SBATCH --partition=production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r2_xsolver_%A.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r2_xsolver_%A.err

# Audit run R2 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 0):
# extend the mesolve vs mcsolve cross-solver benchmark from 2 cells to 3 cells
# (zero and intermediate disorder, both sectors, both loss profiles), so that
# audit item A-P1-12 is demonstrated at three cells instead of assumed.
# Written to a separate JSON so the canonical 2-cell benchmark is preserved.
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

exec "$PY" "$ROOT/scripts/run_cross_solver_benchmark.py" \
  --n-cells 3 --n-max 2 --ntraj 800 \
  --disorders 0.0 0.5 \
  --output "$ROOT/data/cross_solver_benchmark_3cells.json"
