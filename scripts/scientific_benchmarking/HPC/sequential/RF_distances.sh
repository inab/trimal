#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=6
#SBATCH --qos=debug
#SBATCH --output=AA_Prank_RF_computations.txt

module load python
residue_type="AA"
./trimal/scripts/compute_RF_distances.sh $residue_type
