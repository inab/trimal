#!/bin/bash

dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_DNA.txt"

for problem in $dataset/DNA/*/problem*
do
	if [[ $problem != *"problem"*"_err" ]]; then
		problem_name=$(basename $problem | awk -F '.' '{print $1}')
		residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
		if grep -Fxq ${residue_taxon_problem} ${problems_to_ignore}; then
			:
		else
			for file in $problem/*.{ClustalW2,Mafft,Prank,T-Coffee}.fa
			do
				echo -n "./../trimal/scripts/get_statistics.sh "../$file""
				tool=$(basename $file | awk -F '.' '{print $2}')
				for filtered_file in $problem/*.${tool}.*.fa
				do
					if [[ $filtered_file != *"automated2"* ]]; then
						echo -n " && ./../trimal/scripts/get_statistics.sh "../$filtered_file""
					fi
				done
				echo
			done
		fi
	fi
done
