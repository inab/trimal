#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=24
#SBATCH --qos=bsc_ls
#SBATCH --output=output.txt

#module load iq-tree

FILE=file.txt
/apps/GREASY/2.2/INTEL/IMPI/bin/greasy $FILE
