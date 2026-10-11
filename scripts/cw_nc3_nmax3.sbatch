#!/bin/bash
#SBATCH --job-name=cw_nc3_nmax3
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/cw_nc3_nmax3_%j.log
#SBATCH --error=logs/cw_nc3_nmax3_%j.err

# Collective-mode Wigner pilot — truncation-stability extension:
# 3 cells / 6 modes / N_max=3 (Fock dim 4^6 = 4096).
# N_max=2 at 6 modes already measured (data/collective_wigner_dynamics_nc3.json):
# collective edge-mode negativity 0/80, identity exact, local contrast ~0.
set -e
source /home/nanaengo/miniforge3/etc/profile.d/conda.sh
conda activate qom

cd /home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608
mkdir -p logs

python scripts/run_collective_wigner.py \
    --n-cells 3 --n-max 3 --n-mc 10 --ntraj 10 \
    --output data/collective_wigner_dynamics_nc3_nmax3.json \
    --analysis-output data/collective_wigner_analysis_nc3_nmax3.json

echo "CW_NC3_NMAX3_DONE rc=$?"
