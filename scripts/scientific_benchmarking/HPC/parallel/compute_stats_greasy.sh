#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --array=1-297
#SBATCH --qos=bsc_ls
#SBATCH --workdir=compute_stats_AA
#SBATCH --output=compute_stats_AA_%A_%a.txt
#SBATCH --error=compute_stats_AA_%A_%a.err
module load python
$(sed $SLURM_ARRAY_TASK_ID"q;d" generate_stats_calls_greasy.txt)
