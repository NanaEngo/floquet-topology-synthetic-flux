#!/bin/bash
#SBATCH --job-name=trackA_chiral_prod
#SBATCH --partition=production
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_chiral_prod_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_chiral_prod_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
[ -d "$ROOT" ] || { echo "workspace not found: $ROOT" >&2; exit 2; }
mkdir -p "$ROOT/logs" "$ROOT/data/chiral_production"

case "${SLURM_ARRAY_TASK_ID}" in
  0) CELLS=2; N_MC=500 ;;
  1) CELLS=3; N_MC=500 ;;
  2) CELLS=4; N_MC=300 ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

exec "$PY" "$ROOT/scripts/run_chiral_floquet_dynamics.py" \
  --n-cells "$CELLS" \
  --n-max 2 \
  --n-mc "$N_MC" \
  --ntraj 10 \
  --wigner-grid 31 \
  --initial-state-mode edge_site \
  --kerr 0.5 \
  --gamma-edge 1e-4 \
  --gamma-bulk 0.2 \
  --t-stop 6.283185307179586 \
  --n-time 9 \
  --seed 20260814 \
  --output "$ROOT/data/chiral_production/chiral_floquet_${CELLS}cells_NMC${N_MC}.json"
