#!/bin/bash
#SBATCH --job-name=trackA_n6_rebuild
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
ROOT="${SLURM_SUBMIT_DIR}"
PY="/home/nanaengo/miniforge3/envs/qom/bin/python"
export PYTHONPATH="/home/nanaengo/FL_QOM/src:${PYTHONPATH:-}"
mkdir -p "${ROOT}/logs" "${ROOT}/data"
exec > >(tee "${ROOT}/logs/trackA_n6_${SLURM_JOB_ID}.log") 2> >(tee "${ROOT}/logs/trackA_n6_${SLURM_JOB_ID}.err" >&2)

"${PY}" -u "${ROOT}/scripts/run_trackA_n6_independent.py" \
  --n-mc 500 \
  --ntraj 10 \
  --clean-ntraj 50 \
  --seed 42 \
  --output "${ROOT}/data/trackA_n6_independent.json"
