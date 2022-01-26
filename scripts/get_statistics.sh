#!/bin/bash


scripts="/home/nicolas/Desktop/BSC/trimal/scripts/"
dataset="/home/nicolas/Desktop/BSC/trimal/dataset/"
test_working_files="/home/nicolas/Desktop/BSC/test_working_files"
cd $test_working_files
> blocks_outputs.txt
> number_sequences.txt
> number_columns.txt
> gap_stats.txt
for file in $dataset*.fasta;
do
    echo $file
    trimal -in $file -sgc > temporal_out.txt
    if [ -s temporal_out.txt ]; then
        python3 $scripts/set_manual_boundaries.py -i temporal_out.txt --inner_blocks --total_blocks --min_gapscore_allowed 1 >> blocks_outputs.txt
        grep ">" $file | wc -l >> number_sequences.txt
        trimal -in $file -sgc | tail -n 1 | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' >> number_columns.txt
        trimal -in $file -sgt > temporal_out.txt
        gap_score=$(grep -E "[[:digit:]]" temporal_out.txt -m 1 -o | tail -n 1)
        if [ $gap_score == 1 ]; then
            grep -E "[[:digit:]]{1,}" temporal_out.txt -m 1 -o | head -n 1 >> gap_stats.txt
        else
            echo 0 >> gap_stats.txt
        fi
    fi
done