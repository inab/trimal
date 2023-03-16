#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=3
#SBATCH --qos=bsc_ls
#SBATCH --output=tech_bench_greasy_missing.txt
#SBATCH --error=tech_bench_greasy_missing.err
module load python

FILE=tech_bench_array_calls_folder/tech_bench_calls_hq
/apps/GREASY/2.2/INTEL/IMPI/bin/greasy $FILE
