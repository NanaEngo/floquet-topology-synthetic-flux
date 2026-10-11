#!/bin/bash
#SBATCH --job-name=trackA_match_sens
#SBATCH --partition=production
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_match_sens_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_match_sens_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
export PYTHONPATH="/home/nanaengo/FL_QOM/src"
[ -d "$ROOT" ] || { echo "workspace not found: $ROOT" >&2; exit 2; }
mkdir -p "$ROOT/logs" "$ROOT/data/sensitivity"

case "${SLURM_ARRAY_TASK_ID}" in
  0) N=4 ;;
  1) N=6 ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

exec "$PY" "$ROOT/scripts/run_matched_dynamics.py" \
  --n-sites "$N" \
  --n-max 2 \
  --n-mc 10 \
  --ntraj 10 \
  --wigner-grid 31 \
  --gamma-edge 1e-4 \
  --gamma-bulk 0.2 \
  --t-stop 5.0 \
  --n-time 11 \
  --seed 20260814 \
  --output "$ROOT/data/sensitivity/matched_dynamics_N${N}_sensitivity.json"
