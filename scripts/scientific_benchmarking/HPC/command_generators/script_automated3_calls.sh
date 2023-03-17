#!/bin/bash

dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_AA.txt"
scripts="trimal/scripts"

for problem in $dataset/AA/*/problem*
do
	problem_name=$(basename $problem | awk -F '.' '{print $1}')
	residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
	if grep -Fxq ${residue_taxon_problem} ${problems_to_ignore}; then
		:
	else
		for file in $problem/${problem_name}.{ClustalW2,Mafft,Prank,T-Coffee}.fa
		do
			echo "./${scripts}/generate_trees_automated3.sh $file"
		done
	fi
done
