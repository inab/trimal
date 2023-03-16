#!/bin/bash

start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"

line_number=1
for residue_type in $dataset/AA*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem*
            do
                residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
                residue_taxon_problem_to_ignore=$(sed "${line_number}q;d" all_gaps_indets_alignments_unique.txt)
                if [[ $residue_taxon_problem != $residue_taxon_problem_to_ignore ]]; then
                    for alignment in $problem/*.fa
                    do
                        iqtree -nt 2 -quiet -mem 4G -cmin 4 -cmax 10 -s $alignment -bb 1000 -mset WAG,LG,JTT &
                        echo $alignment
                    done
                else
                    ((line_number++))
                fi

            done
		done
    fi
done

wait

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc