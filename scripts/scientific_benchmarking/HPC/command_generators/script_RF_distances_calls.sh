#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --qos=debug
#SBATCH --output=RF_distances_calls_AA_automated3.txt

dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_AA.txt"

for problem in $dataset/AA/**/problem*
do
	problem_name=$(basename $problem | awk -F '_' '{print $1}')
	residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
	if grep -Fxq ${residue_taxon_problem} ${problems_to_ignore}; then
		:
	else
		if [[ $problem != *"problem"*"_err" ]]; then
			for file in $problem/*automated3_*.fa
			do
				echo "./trimal/scripts/compute_RF_distances.sh $file"
			done
		fi
	fi
done
