#!/bin/bash
# Submit from the project root: sbatch script/simulation/submit_simulation.sh
#SBATCH --job-name=susie_sim
#SBATCH --output="simulation results/slurm_%A_%a.out"
#SBATCH --error="simulation results/slurm_%A_%a.err"
#SBATCH --time=23:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-200%20

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script from the susie_mix project root}"
module load R/4.2.0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript --vanilla script/simulation/run_simulation.R \
  --shard-id "${SLURM_ARRAY_TASK_ID}" --n-shards 200 "$@"
