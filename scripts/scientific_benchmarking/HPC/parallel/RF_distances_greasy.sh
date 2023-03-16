#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --array=1-300
#SBATCH --qos=bsc_ls
#SBATCH --output=RF_distances_arrays_AA_automated3/RF_distances_arrays_AA_automated3_%A_%a.txt
#SBATCH --error=RF_distances_arrays_AA_automated3/RF_distances_arrays_AA_automated3_%A_%a.err
module load python
echo "Array task $SLURM_ARRAY_TASK_ID"
$(sed $SLURM_ARRAY_TASK_ID"q;d" RF_distances_list_greasy.txt)
