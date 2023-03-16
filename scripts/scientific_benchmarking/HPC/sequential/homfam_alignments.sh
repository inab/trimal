#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=debug
#SBATCH --output=homfam_alignments.txt
module load python
./trimal/scripts/filter_homfam_alignments.sh
