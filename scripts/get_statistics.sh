#!/bin/bash


scripts="/home/nicolas/Desktop/BSC/trimal/scripts/"
dataset="/home/nicolas/Desktop/BSC/trimal/dataset/"
> blocks_outputs.txt
for file in $dataset*.fasta;
do
    echo $file
    trimal -in $file -sgc > temporal_out.txt
    [ -s temporal_out.txt ] && python3 $scripts/set_manual_boundaries.py -i temporal_out.txt --inner_blocks --total_blocks --min_gapscore_allowed 1 >> blocks_outputs.txt
done
#python3 trimal/scripts/set_manual_boundaries.py -i Outputs/example.015.AA.txt --inner_blocks --total_blocks