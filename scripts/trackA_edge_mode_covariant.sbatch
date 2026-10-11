#!/bin/bash
#SBATCH --job-name=trackA_edge_cov
#SBATCH --partition=production
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_edge_cov_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_edge_cov_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
[ -d "$ROOT" ] || { echo "workspace not found: $ROOT" >&2; exit 2; }
mkdir -p "$ROOT/logs" "$ROOT/data/edge_mode_covariant"

case "${SLURM_ARRAY_TASK_ID}" in
  0) PREP=common_edge_site; KERR=0.0; TAG=common_passive ;;
  1) PREP=common_edge_site; KERR=0.5; TAG=common_kerr ;;
  2) PREP=topological_edge_mode; KERR=0.0; TAG=mode_passive ;;
  3) PREP=topological_edge_mode; KERR=0.5; TAG=mode_kerr ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

exec "$PY" "$ROOT/scripts/run_chiral_edge_mode_pilot.py" \
  --n-cells 2 \
  --n-max 2 \
  --n-mc 10 \
  --ntraj 5 \
  --wigner-grid 31 \
  --initial-state-mode "$PREP" \
  --edge-observable-frame micromotion \
  --edge-operator-convention standard \
  --kerr "$KERR" \
  --gamma-edge 1e-4 \
  --gamma-bulk 0.2 \
  --t-stop 50.26548245743669 \
  --n-time 17 \
  --seed 20260819 \
  --output "$ROOT/data/edge_mode_covariant/edge_mode_${TAG}.json"
