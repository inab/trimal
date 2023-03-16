#!/bin/bash

scripts="trimal/scripts"
file=$1
method=$2
min_percentage_columns=$3
min_columns=$4


task(){
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    folder_path=$(dirname $file)
	file_basename=$(basename $file .fa)
	output_file_prefix="$folder_path/${file_basename}"

	echo "Calculating trimmed alignment for ${file_basename}"

	# generate trimmed alignments
	(python3 $scripts/automated2_algorithm.py -i $file -a $method --min_percentage_columns $min_percentage_columns \
		--min_columns $min_columns > ${output_file_prefix}.trimAl_automated2_${method}_${min_percentage_columns}_${min_columns}.fa)

	echo "Calculating tree for ${file_basename}.trimAl_automated2.fa"
	# generate new trees
	iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s ${output_file_prefix}.trimAl_automated2_${method}_${min_percentage_columns}_${min_columns}.fa -bb 1000 -mset GTR
}


task "$file"
