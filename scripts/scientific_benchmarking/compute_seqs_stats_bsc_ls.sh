#!/bin/bash

start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"

> all_gaps_indets_alignments.txt

calculate_all_gaps_indets_alignments() {
    residue_type=$1
    alignment=$2
    residue_type_argument=$(echo $residue_type | awk -F '/' '{print $NF}')
    all_gaps_indets=$(python3 trimal/scripts/get_alignment_statistics.py -i $alignment -t $residue_type_argument --all_gaps_indets)
    if [[ $all_gaps_indets == "True" ]]; then
        echo $alignment >> all_gaps_indets_alignments.txt
    fi
}

counter=0
for residue_type in $dataset/*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem*
            do
                for alignment in $problem/*.fa
                do
                    calculate_all_gaps_indets_alignments $residue_type $alignment &
		    echo $alignment
		    ((counter++))
		    if [ $counter -eq 48 ]; then
		        wait
			counter=0
		    fi
                done
            done
		done
    fi
done

wait

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc
