#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=6
#SBATCH --qos=bsc_ls
#SBATCH --output=AA_Prank_new_trees_computations.txt

module load python
residue_type="AA"
./trimal/scripts/get_statistics.sh $residue_type
