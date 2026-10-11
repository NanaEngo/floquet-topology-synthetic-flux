#!/bin/bash
#SBATCH --job-name=trackA_passive_1p
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_passive_1p_%j.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_passive_1p_%j.err

set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
cd "$ROOT"
"$PY" scripts/run_passive_single_particle_benchmark.py \
  --sizes 2 3 4 6 \
  --profiles edge_low_loss uniform_loss \
  --n-time 257 \
  --output data/passive_single_particle_benchmark.json
