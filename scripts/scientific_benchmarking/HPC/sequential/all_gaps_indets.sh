#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=debug
#SBATCH --output=all_gaps_indets.txt
module load python
./trimal/scripts/compute_seqs_stats.sh
