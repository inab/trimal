#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=bsc_ls
#SBATCH --output=homfam_SoP_greasy_seqs.txt
module load python

echo "file,clean_num_seqs,trimmed_num_seqs_gappyout,trimmed_num_seqs_strictplus" > homfam_SoP_stats_seqs.txt

FILE=homfam_SoP_calls.txt
/apps/GREASY/2.2/INTEL/IMPI/bin/greasy $FILE >> homfam_SoP_stats_seqs.txt
