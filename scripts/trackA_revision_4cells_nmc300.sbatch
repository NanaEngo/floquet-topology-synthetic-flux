#!/bin/bash
#SBATCH --job-name=trackA_rev_4c_nmc300
#SBATCH --partition=production
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_rev_4c_nmc300_%A.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_rev_4c_nmc300_%A.err

set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

exec "$PY" "$ROOT/scripts/run_revision_scans.py" gamma \
  --n-cells 4 --n-mc 300 \
  --ratios 5e-4 \
  --disorder-levels 0.0 0.5 1.0 1.5 \
  --seed 20260814 \
  --output "$ROOT/data/chiral_production/chiral_floquet_4cells_NMC300.json"