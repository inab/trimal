#!/bin/bash

start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"

for residue_type in $dataset/DNA*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem*
            do
                for alignment in $problem/*Prank*.fa
                do
                    echo $alignment
                    iqtree -nt 2 -quiet -mem 4G -cmin 4 -cmax 10 -s $alignment -bb 1000 -mset GTR &
                done
            done
		done
    fi
done

wait

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc
