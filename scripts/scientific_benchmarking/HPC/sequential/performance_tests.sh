#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=debug
#SBATCH --output=performance_tests.txt
./trimal/scripts/technical_performance_tests.sh
