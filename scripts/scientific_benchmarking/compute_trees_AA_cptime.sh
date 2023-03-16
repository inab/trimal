#!/bin/bash
res_type=$1
problem_nums=$2
taxon_name=$3
start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_$res_type.txt"
grep "${taxon_name}_problem$problem_nums" $problems_to_ignore > "all_gaps_indets_alignments_unique_${res_type}_${taxon_name}_${problem_nums}.txt"
problems_to_ignore="all_gaps_indets_alignments_unique_${res_type}_${taxon_name}_${problem_nums}.txt"

counter=0
line_number=1
for residue_type in $dataset/$res_type*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/$taxon_name*
        do
            for problem in $taxon/problem$problem_nums*
            do
                residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
                residue_taxon_problem_to_ignore=$(sed "${line_number}q;d" $problems_to_ignore)
                if [[ $residue_taxon_problem != $residue_taxon_problem_to_ignore ]]; then
                    for alignment in $problem/*.fa
                    do
                        if [[ $alignment == *"automated2"* ]]; then
                            iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $alignment -bb 1000 -mset WAG,LG,JTT &
                            echo $alignment
                        fi
                    done
                else
                    ((line_number++))
                fi
                ((counter++))
                if [ $counter -eq 2 ]; then
                    wait
                    counter=0
                fi
            done
		done
    fi
done

wait

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc
