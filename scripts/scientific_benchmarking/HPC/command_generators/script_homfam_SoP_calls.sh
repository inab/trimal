#!/bin/bash

dataset_alignments="cedric_datasets/alignments"
> homfam_SoP_calls.txt
for file in $dataset_alignments/*.aln
do
	echo "./trimal/scripts/filter_homfam_alignments.sh $file" >> homfam_SoP_calls.txt
done
