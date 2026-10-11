#!/bin/bash
#SBATCH --job-name=trackA_r8_witness
#SBATCH --partition=production
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r8_witness_%A_%a.out
#SBATCH --error=/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/logs/trackA_r8_witness_%A_%a.err

# Audit run R8 (AUDIT_RESEARCH_REFINEMENTS_20261001.md, Tier 1): bit-exact
# rerun of the production design (seed 20260814, N_MC=500, all four disorder
# levels, both profiles, both sectors) at 2/3/4 cells, dumping the terminal
# edge/centre reduced blocks.  analyze_r8_witness.py fails closed unless the
# rerun negativities reproduce the archived production records exactly, then
# computes the Chabaud Fock witnesses and the relative-entropy
# non-Gaussianity on the identical trajectory ensemble.
set -euo pipefail
ROOT="/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
mkdir -p "$ROOT/logs" "$ROOT/data/r8_witness"

case "${SLURM_ARRAY_TASK_ID}" in
  0) CELLS=2 ;;
  1) CELLS=3 ;;
  2) CELLS=4 ;;
  *) echo "unexpected array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$PY" "$ROOT/scripts/run_r3_r8_records.py" \
  --n-cells "$CELLS" \
  --n-mc 500 \
  --dump-blocks \
  --n-workers 16 \
  --seed 20260814 \
  --output "$ROOT/data/r8_witness/r8_witness_${CELLS}cells.json"
