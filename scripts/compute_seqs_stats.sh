#!/bin/bash

start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
action=$1

calculate_all_gaps_indets_alignments() {
    residue_type=$1
    alignment=$2
    residue_type_argument=$(echo $residue_type | awk -F '/' '{print $NF}')
    all_gaps_indets=$(python3 trimal/scripts/get_alignment_statistics.py -i $alignment -t $residue_type_argument --all_gaps_indets)
    if [[ $all_gaps_indets == "True" ]]; then
        echo $alignment >> all_gaps_indets_alignments.txt
    fi
}

if [[ $action == "calculate_all_gaps_indets" ]]; then
    > all_gaps_indets_alignments.txt
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
fi

wait


if [[ $action == "plot_sequences_statistics" ]]; then
    > empty_sequences_statistics.csv
    variable_names=()
    variable_values=()
    residue_types=("AA" "DNA")
    taxa=("Bacteria" "Eukaryotes" "Fungi")
    MSA_tools=("ClustalW2" "ClustalW" "Mafft" "Prank" "T-Coffee")
    MSA_filter_tools=("Aliscore" "BMGE62" "BMGE100" "Gblocks" "Noisy" "trimAl" "Zorro" "Guidance")
    for residue_type in ${residue_types[@]}
    do
        echo $residue_type
        variable_names+=(${residue_type}_empty_seqs)
        variable_values+=($(grep -c "$residue_type" all_gaps_indets_alignments.txt))
        for taxon in ${taxa[@]}
        do
            echo "${residue_type}_${taxon}"
            variable_names+=(${residue_type}_${taxon}_empty_seqs)
            variable_values+=($(grep -c "${residue_type}/${taxon}" all_gaps_indets_alignments.txt))
            for MSA_tool in ${MSA_tools[@]}
            do
                echo "${residue_type}_${taxon}_${MSA_tool}"
                variable_names+=(${residue_type}_${taxon}_${MSA_tool}_empty_seqs)
                variable_values+=($(grep -c "${residue_type}/${taxon}/.*${MSA_tool}.*" all_gaps_indets_alignments.txt))
                for MSA_filter_tool in ${MSA_filter_tools[@]}
                do
                    echo "${residue_type}_${taxon}_${MSA_tool}"
                    variable_names+=(${residue_type}_${taxon}_${MSA_tool}_empty_seqs)
                    variable_values+=($(grep -c "${residue_type}/${taxon}/.*${MSA_tool}.*" all_gaps_indets_alignments.txt))
                done
            done
        done
    done
    printf -v joined '%s,' "${variable_names[@]}"
    echo "${joined%,}" >> empty_sequences_statistics.csv
    printf -v joined '%s,' "${variable_values[@]}"
    echo "${joined%,}" >> empty_sequences_statistics.csv
fi

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc
