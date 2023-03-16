#!/bin/bash

start_time=$(date +%s.%N)
python_script="trimal/scripts/encode_columns.py"
dataset="trimal/dataset"

for alignment in $dataset/*.fasta
do
	filename=$(basename $alignment .fasta)
	(python3 $python_script -i $alignment > "encoded_alignments/$filename.txt" && echo $filename) &
done
wait

end_time=$(date +%s.%N)
echo "$end_time - $start_time" | bc
