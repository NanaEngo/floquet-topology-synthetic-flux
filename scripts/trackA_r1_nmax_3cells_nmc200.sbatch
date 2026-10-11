#!/bin/bash
#SBATCH --job-name=trackA_r1_nmax_3c
#SBATCH --partition=production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r1_nmax_3c_%A.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r1_nmax_3c_%A.err

# Audit run R1 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 0):
# 3-cell higher-cutoff sign-consistency panel raised from N_MC=100 (N_max=3)
# and N_MC=50 (N_max=4) to N_MC=200 per cutoff, closing gate G3.
# Fresh output paths are used on purpose: run_revision_scans.py resumes by
# skipping pre-existing fragments, so the historical 100/50-realization
# fragment directories must not be reused.
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

for NMAX in 3 4; do
  "$PY" "$ROOT/scripts/run_revision_scans.py" nmax \
    --n-cells 3 --n-max "$NMAX" --n-mc 200 \
    --disorder-levels 0.0 0.35 \
    --seed 20260817 \
    --output "$ROOT/data/revision_scans/nmax_${NMAX}_3cells_nmc200.json"
done
