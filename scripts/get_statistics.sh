#!/bin/bash

file=$1
residue_type="DNA"
trimal_local="../trimal/bin/trimal"
scripts="../trimal/scripts"

task() {
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
	folder_path=$(dirname $file)
	if [[ $file == *".fa" ]]; then
		tree_file="$file.treefile"
	fi
	if [[ $filename != *".seqs."* ]]; then
		$trimal_local -in $file -sgc > temp_out_$filename.txt
		(python3 $scripts/set_manual_boundaries.py -i temp_out_$filename.txt --inner_blocks --total_blocks --block_coordinates \
			--min_gapscore_allowed 1 --max_blocks 10000 > temp_number_blocks_$filename.txt)
		number_blocks=$(grep "## Blocks" temp_number_blocks_$filename.txt | awk '{print $3}')
		left_block_column=$(grep "## Left column" temp_number_blocks_$filename.txt | awk '{print $4}')
		right_block_column=$(grep "## Right column" temp_number_blocks_$filename.txt | awk '{print $4}')
		rm temp_number_blocks_$filename.txt
		number_columns=$(tail -n 1 temp_out_$filename.txt | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}')
		rm temp_out_$filename.txt
		$trimal_local -in $file -sgt > temp_out_sgt_$filename.txt
		avg_gaps=$(grep "## AverageGaps" temp_out_sgt_$filename.txt | awk '{print $3}')
		rm temp_out_sgt_$filename.txt
		avg_seq_identity=$($trimal_local -in $file -sident | grep "## AverageIdentity" | awk '{print $3}')
		if [ -s $tree_file ]; then
			ref_tree=$(find $folder_path -name *.reference.nwk)
			RF_distance=$(ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename.txt &&
			awk 'FNR == 3 {print $9}' temp_RF_distance_$filename.txt || 
			echo -1)
		fi
	fi 

    if [ -s $file.treefile ]; then
        rm temp_RF_distance_$filename.txt
    fi

	num_seq=$(grep -c ">" $file)
	taxon=$(echo $filename | awk -F '_' '{print $2}')
	problem_num=$(echo $filename | awk -F '.' '{print $1}' | grep -o '[0-9]\+')
	MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
	MSA_filter_tool=""
	if [[ $MSA_tool == Guidance* ]]; then
		MSA_tool=$(echo $MSA_tool | awk -F 'Guidance' '{print $2}')
		MSA_filter_tool='Guidance'
	else
		MSA_filter_tool=$(echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}')
	fi

	residue_taxon_problem=$(echo $filename | awk -F '.' '{print $1}')

	if [[ $MSA_filter_tool != "None" ]]; then
		if [[ $MSA_filter_tool == "Guidance" ]]; then
			if [[ $MSA_tool == "ClustalW" ]]; then
				avg_gaps_MSA_tool_file="avg_gaps_$residue_taxon_problem."$MSA_tool"2.fa.txt"
				avg_seq_identity_MSA_tool_file="avg_seq_identity_$residue_taxon_problem."$MSA_tool"2.fa.txt"
				number_columns_MSA_tool_file="number_columns_$residue_taxon_problem."$MSA_tool"2.fa.txt"
				RF_distance_MSA_tool_file="RF_distance_$residue_taxon_problem."$MSA_tool"2.fa.txt"
			else
				avg_gaps_MSA_tool_file="avg_gaps_$residue_taxon_problem.$MSA_tool.fa.txt"
				avg_seq_identity_tool_file="avg_seq_identity_$residue_taxon_problem.$MSA_tool.fa.txt"
				number_columns_MSA_tool_file="number_columns_$residue_taxon_problem.$MSA_tool.fa.txt"
				RF_distance_MSA_tool_file="RF_distance_$residue_taxon_problem.$MSA_tool.fa.txt"
			fi
		else
			file_suffix=$(echo $filename | awk -F '.' '{print $(NF-3)"."$(NF-2)"."$NF}')
			avg_gaps_MSA_tool_file="avg_gaps_$file_suffix.txt"
			avg_seq_identity_MSA_tool_file="avg_seq_identity_$file_suffix.txt"
			number_columns_MSA_tool_file="number_columns_$file_suffix.txt"
			RF_distance_MSA_tool_file="RF_distance_$file_suffix.txt"
		fi
		avg_gaps_MSA_tool=$(cat $avg_gaps_MSA_tool_file)
		avg_gaps_diff=$(echo "$avg_gaps_MSA_tool - $avg_gaps" | bc)
		avg_seq_identity_MSA_tool=$(cat $avg_seq_identity_MSA_tool_file)
		avg_seq_identity_diff=$(echo "$avg_seq_identity_MSA_tool - $avg_seq_identity" | bc)
		number_columns_tool=$(cat $number_columns_MSA_tool_file)
		removed_columns=$(echo "$number_columns_tool - $number_columns" | bc)
		percent_conserved_columns=$(echo "$number_columns / $number_columns_tool" | bc -l)
		RF_distance_tool=$(cat $RF_distance_MSA_tool_file)
		RF_distance_diff=$(echo "$RF_distance_tool - $RF_distance" | bc)
		avg_gaps_diff_weighted=$(echo "$avg_gaps_diff * $percent_conserved_columns" | bc -l)
		avg_seq_identity_diff_weighted=$(echo "$avg_seq_identity_diff * $percent_conserved_columns" | bc -l)
	else
		(echo $avg_gaps || echo -1 ) > avg_gaps_$filename.txt
		(echo $avg_seq_identity || echo -1 ) > avg_seq_identity_$filename.txt
		(echo $RF_distance || echo -1 ) > RF_distance_$filename.txt
		(echo $number_columns || echo -1 ) > number_columns_$filename.txt
	fi

	echo "$residue_type,$taxon,$problem_num,$MSA_tool,$MSA_filter_tool,$number_columns_tool,$num_seq,$number_columns,$number_blocks,$left_block_column,$right_block_column,$avg_gaps,$avg_seq_identity,$RF_distance,$removed_columns,$percent_conserved_columns,$avg_gaps_diff_weighted,$avg_seq_identity_diff_weighted,$RF_distance_diff"
}


task "$file"