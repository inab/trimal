#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=debug
#SBATCH --output=RF_distances_greasy.txt
module load python

echo "taxon,problem_num,number_sequences,number_columns,MSA_tool,MSA_filter_tool,RF_distance,RF_distance_diff" > RF_distances_all_tools.txt
FILE=RF_distances_calls.txt
/apps/GREASY/2.2/INTEL/IMPI/bin/greasy $FILE >> RF_distances_all_tools.txt
sed -i "2d"  RF_distances_all_tools.txt
