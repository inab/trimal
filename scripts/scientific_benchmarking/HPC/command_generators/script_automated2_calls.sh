#!/bin/bash

dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_DNA.txt"
scripts="trimal/scripts"
min_percentage_columns_list=(30 35 40 45)
min_columns_list=(100 150 200 250)
methods=("gaps" "similarity" "combined")

for problem in $dataset/DNA/*/problem*
do
	problem_name=$(basename $problem | awk -F '.' '{print $1}')
	residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
	if grep -Fxq ${residue_taxon_problem} ${problems_to_ignore}; then
		:
	else
		for file in $problem/${problem_name}.{ClustalW2,Mafft,Prank,T-Coffee}.fa
		do
			for min_percentage_columns in ${min_percentage_columns_list[*]}; do
				for min_columns in ${min_columns_list[*]}; do
					for method in ${methods[*]}; do
						echo "./${scripts}/generate_trees_automated2.sh $file $method $min_percentage_columns $min_columns"
					done
				done
			done
		done
	fi
done
