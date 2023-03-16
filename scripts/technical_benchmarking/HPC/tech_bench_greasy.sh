#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=2
#SBATCH --qos=bsc_ls
#SBATCH --array=1-300
#SBATCH --workdir=tech_bench_files
#SBATCH --output=tech_bench_greasy_chunk1_from_5GB_to_20GB_%A_%a.txt
#SBATCH --error=tech_bench_greasy_chunk1_from_5GB_to_20GB_%A_%a.err

FILE=tech_bench_calls_greasy_chunk1_from_5GB_to_20GB.txt
$(sed $SLURM_ARRAY_TASK_ID"q;d" $FILE);
