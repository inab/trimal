#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=24
#SBATCH --array=1-300
#SBATCH --qos=bsc_ls
#SBATCH --output=compute_trees_DNA/compute_trees_DNA_%A_%a.txt
#SBATCH --error=compute_trees_DNA/compute_trees_DNA_%A_%a.err
module load iq-tree
$(sed $SLURM_ARRAY_TASK_ID"q;d" generate_trees_calls_greasy.txt)
