#!/bin/bash
#SBATCH --job-name=trackA_rev_nmax4_3c
#SBATCH --partition=production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_rev_nmax4_3c_%A.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_rev_nmax4_3c_%A.err

set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

exec "$PY" "$ROOT/scripts/run_revision_scans.py" nmax \
  --n-cells 3 --n-max 4 --n-mc 50 \
  --disorder-levels 0.0 0.35 \
  --seed 20260817 \
  --output "$ROOT/data/revision_scans/nmax_4_3cells.json"