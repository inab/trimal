#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=bsc_ls
#SBATCH --output=all_gaps_indets_bsc_ls.txt
module load python
./trimal/scripts/compute_seqs_stats_bsc_ls.sh
