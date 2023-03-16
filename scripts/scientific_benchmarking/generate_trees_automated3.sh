#!/bin/bash

scripts="trimal/scripts"
trimal_local="trimal/bin/trimal"
file=$1


task(){
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    folder_path=$(dirname $file)
	file_basename=$(basename $file .fa)
	output_file_prefix="$folder_path/${file_basename}"

	echo "Calculating trimmed alignment for ${file_basename}"

	gappyout_cols=$(awk '/^>/{l=0; next}{l+=length($0)}END{print l}' ${output_file_prefix}.trimAl_gappyout.fa)
	trimAl_automated2_gaps_30_100_cols=$(awk '/^>/{l=0; next}{l+=length($0)}END{print l}' ${output_file_prefix}.trimAl_automated2_gaps_30_100.fa)
	if (( $(echo "$gappyout_cols > $trimAl_automated2_gaps_30_100_cols" | bc -l) )) ; then
		cat ${output_file_prefix}.trimAl_gappyout.fa > ${output_file_prefix}.trimAl_automated3_1.fa
		echo "gappyout $gappyout_cols bigger than automated2_gaps_30_100 $trimAl_automated2_gaps_30_100_cols"
	else
		cat ${output_file_prefix}.trimAl_automated2_gaps_30_100.fa > ${output_file_prefix}.trimAl_automated3_1.fa
		echo "automated2_gaps_30_100 $trimAl_automated2_gaps_30_250_cols bigger than gappyout $gappyout_cols"
	fi

	trimAl_automated2_gaps_30_250_cols=$(awk '/^>/{l=0; next}{l+=length($0)}END{print l}' ${output_file_prefix}.trimAl_automated2_gaps_30_250.fa)
	if (( $(echo "$gappyout_cols > $trimAl_automated2_gaps_30_250_cols" | bc -l) )) ; then
		cat ${output_file_prefix}.trimAl_gappyout.fa > ${output_file_prefix}.trimAl_automated3_2.fa
		echo "gappyout $gappyout_cols bigger than automated2_gaps_30_250 $trimAl_automated2_gaps_30_250_cols"
	else
		cat ${output_file_prefix}.trimAl_automated2_gaps_30_250.fa > ${output_file_prefix}.trimAl_automated3_2.fa
		echo "automated2_gaps_30_250 $trimAl_automated2_gaps_30_250_cols bigger than gappyout $gappyout_cols"
	fi


	# generate new trees
	iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s ${output_file_prefix}.trimAl_automated3_1.fa -bb 1000 -mset WAG,LG,JTT
	iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s ${output_file_prefix}.trimAl_automated3_2.fa -bb 1000 -mset WAG,LG,JTT
}


task "$file"
