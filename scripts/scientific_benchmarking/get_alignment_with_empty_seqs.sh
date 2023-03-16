#!/bin/bash

start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"

> all_gaps_indets_alignments.txt

for residue_type in $dataset/*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem*
            do
                for alignment in $problem/*Prank*.fa
                do
                    (residue_type_argument=$(echo $residue_type | awk -F '/' '{print $NF}')
                    if [[ $(python3 trimal/scripts/get_alignment_statistics.py -i $alignment -t $residue_type_argument --all_gaps_indets) == "True" ]]; then
                        echo $alignment >> all_gaps_indets_alignments.txt
                    fi
                    echo $alignment) &
                done
            done
		done
    fi
done

wait

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc