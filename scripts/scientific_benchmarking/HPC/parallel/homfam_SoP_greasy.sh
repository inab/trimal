#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=4
#SBATCH --array=1-277
#SBATCH --qos=bsc_ls
#SBATCH --output=homfam_SoP_folder/homfam_SoP_arrays_%A_%a.txt
#SBATCH --error=homfam_SoP_folder/homfam_SoP_arrays_%A_%a.err
module load python
$(sed $SLURM_ARRAY_TASK_ID"q;d" homfam_SoP_list_greasy.txt );
