#!/bin/bash
#SBATCH --job-name=trackA_nmc2000
#SBATCH --partition=production
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_nmc2000_%j.log
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_nmc2000_%j.err

# Power extension: 8-mode (4 cells) at N_MC=2000 (4x the N_MC=500 production,
# nested samples 0-499 identical by seed construction). Row-level checkpoint
# via run_8mode_nmc500_parallel.py (.partial.jsonl resume on requeue).
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
mkdir -p "$ROOT/logs" "$ROOT/data/chiral_production"

echo "=== 4cells N_MC=2000 start $(date -Is)"
"$PY" "$ROOT/scripts/run_8mode_nmc500_parallel.py" \
  --n-cells 4 \
  --n-max 2 \
  --n-mc 2000 \
  --ntraj 10 \
  --n-workers 16 \
  --wigner-grid 31 \
  --kerr 0.5 \
  --gamma-edge 1e-4 \
  --gamma-bulk 0.2 \
  --t-stop 6.283185307179586 \
  --n-time 9 \
  --seed 20260814 \
  --output "$ROOT/data/chiral_production/chiral_floquet_4cells_NMC2000.json"
echo "NMC2000_8MODE_DONE rc=$?"
