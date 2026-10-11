#!/bin/bash
#SBATCH --job-name=trackA_mapped_pilot
#SBATCH --partition=production
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_mapped_pilot_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_mapped_pilot_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
export PYTHONPATH="/home/nanaengo/FL_QOM/src"
[ -d "$ROOT" ] || { echo "workspace not found: $ROOT" >&2; exit 2; }
mkdir -p "$ROOT/logs" "$ROOT/data/mapped_pilot"

case "${SLURM_ARRAY_TASK_ID}" in
  0) CELLS=2 ;;
  1) CELLS=3 ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

exec "$PY" "$ROOT/scripts/run_mapped_floquet_dynamics.py" \
  --n-cells "$CELLS" \
  --n-max 2 \
  --n-mc 10 \
  --ntraj 10 \
  --wigner-grid 31 \
  --gamma-edge 1e-4 \
  --gamma-bulk 0.2 \
  --t-stop 6.283185307179586 \
  --n-time 9 \
  --seed 20260814 \
  --output "$ROOT/data/mapped_pilot/mapped_floquet_${CELLS}cells.json"
