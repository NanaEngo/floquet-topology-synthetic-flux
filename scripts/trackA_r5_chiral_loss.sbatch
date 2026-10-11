#!/bin/bash
#SBATCH --job-name=trackA_r5_chiral
#SBATCH --partition=production
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r5_chiral_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r5_chiral_%A_%a.err

# Audit run R5 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 1):
# sublattice-localized dissipation tested against the "symmetry-respecting or
# engineered loss" row of Table I, currently marked "Not tested".
#   edge_low_loss             : control (generic on-site loss, [S, a_j] != 0 on odd j)
#   chiral_even_loss          : loss only on the even (A) sublattice -> [S, L_k] = 0
#                               for every collapse operator (strong chiral symmetry)
#   chiral_odd_loss           : loss only on the odd (B) sublattice -> {S, L_k} = 0
#                               (dissipator still invariant, operators anticommuting)
#   chiral_even_loss_matched  : even-sublattice loss rescaled by x2 so that the total
#                               rate budget and the edge/bulk ratio match the control
#                               (closes the "you halved the loss" objection)
# Same matched protocol as production: paired disorder draws, nested
# trajectories, edge_site even-cat preparation, chi=0.5, 31x31 Wigner grid.
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"

CELLS=$(( SLURM_ARRAY_TASK_ID + 2 ))

exec "$PY" "$ROOT/scripts/run_revision_scans.py" chiralloss \
  --n-cells "$CELLS" \
  --n-mc 500 \
  --ntraj 10 \
  --disorder-levels 0.0 0.35 \
  --seed 20260817 \
  --output-dir "$ROOT/data/revision_scans"
