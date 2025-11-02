#!/bin/bash

#SBATCH --mem=32GB
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=./logs/%a_%A.out
#SBATCH --error=./logs/%a_%A.err
#SBATCH --account=pr_181_general

#SBATCH --job-name=FWB_with_CF
#SBATCH --mail-type=END
#SBATCH --mail-user=sd3927@nyu.edu


module purge
module load matlab/2023a

echo "$SLURM_ARRAY_TASK_ID"

matlab -batch "SimulationMain, exit"
