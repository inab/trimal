#!/bin/bash

dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_DNA.txt"

for problem in $dataset/DNA/*/problem*
do
	problem_name=$(basename $problem | awk -F '.' '{print $1}')
	residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
	if grep -Fxq ${residue_taxon_problem} ${problems_to_ignore}; then
		:
	else
		for file in $problem/*.{ClustalW2,Mafft,Prank,T-Coffee}*.fa
		do
			echo "iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $file -bb 1000 -mset GTR"
		done
	fi
done
