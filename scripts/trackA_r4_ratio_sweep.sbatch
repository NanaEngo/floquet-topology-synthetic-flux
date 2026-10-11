#!/bin/bash
#SBATCH --job-name=trackA_r4_sweep
#SBATCH --partition=production
#SBATCH --array=0-37
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r4_sweep_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r4_sweep_%A_%a.err

# Audit run R4 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 1):
# continuous intra/inter sweep at fixed intra+inter = 1.3, i.e. the production
# pair (0.1, 1.2) and (1.2, 0.1) are the endpoints of the path.  The grid is
# refined around the three bulk critical points of the path:
#   pi-gap closes at intra = 0.15 and intra = 1.15 (|ja-jb| = 1 at k = pi),
#   zero-gap closes at intra = 0.65 (ja = jb).
# One drive per point; the observable is the edge-minus-centre negativity
# contrast of the same run, so no mirror-symmetric (antisymmetric by
# construction) pairing is used.  Disorder draws are paired across intra
# values (draw_seed does not depend on the sweep coordinate).
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

INTRAS=(0.1 0.13 0.15 0.17 0.2 0.4 0.5 0.6 0.625 0.65 0.675 0.7 0.8 1.0 1.1 1.13 1.15 1.17 1.2)
N=${#INTRAS[@]}

CELLS=$(( SLURM_ARRAY_TASK_ID / N + 2 ))
INTRA=${INTRAS[$(( SLURM_ARRAY_TASK_ID % N ))]}

exec "$PY" "$ROOT/scripts/run_revision_scans.py" sweep \
  --n-cells "$CELLS" \
  --intra "$INTRA" \
  --total 1.3 \
  --n-mc 500 \
  --ntraj 10 \
  --disorder-levels 0.0 0.35 \
  --seed 20260817 \
  --output-dir "$ROOT/data/revision_scans"
