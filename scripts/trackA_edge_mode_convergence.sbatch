#!/bin/bash
#SBATCH --job-name=trackA_edge_conv
#SBATCH --partition=production
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_edge_conv_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_edge_conv_%A_%a.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
[ -d "$ROOT" ] || { echo "workspace not found: $ROOT" >&2; exit 2; }
mkdir -p "$ROOT/logs" "$ROOT/data/edge_mode_convergence"

case "${SLURM_ARRAY_TASK_ID}" in
  0) NMAX=2; KERR=0.0; TAG=nmax2_passive ;;
  1) NMAX=3; KERR=0.0; TAG=nmax3_passive ;;
  2) NMAX=4; KERR=0.0; TAG=nmax4_passive ;;
  3) NMAX=2; KERR=0.5; TAG=nmax2_kerr ;;
  4) NMAX=3; KERR=0.5; TAG=nmax3_kerr ;;
  5) NMAX=4; KERR=0.5; TAG=nmax4_kerr ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

exec "$PY" "$ROOT/scripts/run_chiral_edge_mode_pilot.py" \
  --n-cells 2 \
  --n-max "$NMAX" \
  --n-mc 3 \
  --ntraj 3 \
  --wigner-grid 31 \
  --initial-state-mode common_edge_site \
  --profiles edge_low_loss \
  --disorder-levels 0.0 \
  --kerr "$KERR" \
  --gamma-edge 1e-4 \
  --gamma-bulk 0.2 \
  --t-stop 50.26548245743669 \
  --n-time 17 \
  --seed 20260818 \
  --output "$ROOT/data/edge_mode_convergence/edge_mode_${TAG}.json"
