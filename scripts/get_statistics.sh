#!/bin/bash

res_type=$1
start_time=$(date +%s.%N)
scripts="../trimal/scripts"
dataset="../dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
alignment_statistics="alignment_statistics_all_tools_$res_type"
trimal_local="../trimal/bin/trimal"
readal_local="../trimal/bin/readal"
problems_to_ignore="../all_gaps_indets_alignments_unique_$res_type.txt"
if [ ! -d "$alignment_statistics" ]; then
    mkdir $alignment_statistics
fi
cd $alignment_statistics

> number_blocks.txt
> left_block_column.txt
> right_block_column.txt
> blocks_diff.txt
> number_sequences.txt
> number_columns.txt
> number_columns_MSA.txt
> max_columns.txt
> min_columns.txt
> removed_columns.txt
> percent_conserved_columns.txt
> avg_gaps.txt
> avg_gaps_diff.txt
> avg_gaps_diff_weighted.txt
> avg_seq_identity.txt
> avg_seq_identity_diff.txt
> avg_seq_identity_diff_weighted.txt
> gappy_columns_50.txt
> gappy_columns_80.txt
> RF_distance.txt
> RF_distance_diff.txt
> residue_type.txt
> taxon.txt
> problem_num.txt
> error_problem.txt
> MSA_tools.txt
> MSA_filter_tools.txt
> log.txt
> filenames_to_write.txt


task(){
    file=$1
    start_file_time=$(date +%s.%N)
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    if [ ! -s temp_out_$filename.txt ]; then
        folder_path=$(dirname $file)
        if [[ $file == *".fa" ]]; then
            tree_file="$file.treefile"
        fi
        if [[ $filename != *".seqs."* ]]; then
            $trimal_local -in $file -sgc > temp_out_$filename.txt
            (python3 $scripts/set_manual_boundaries.py -i temp_out_$filename.txt --inner_blocks --total_blocks --block_coordinates \
			 --min_gapscore_allowed 1 --max_blocks 10000 > temp_number_blocks_$filename.txt && 
			grep "## Blocks" temp_number_blocks_$filename.txt | awk '{print $3}' > number_blocks_$filename.txt &&
			grep "## Left column" temp_number_blocks_$filename.txt | awk '{print $4}' > left_block_column_$filename.txt &&
			grep "## Right column" temp_number_blocks_$filename.txt | awk '{print $4}' > right_block_column_$filename.txt &&
			rm temp_number_blocks_$filename.txt) &
            tail -n 1 temp_out_$filename.txt | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' > number_columns_$filename.txt &
			($trimal_local -in $file -sgt > temp_out_sgt_$filename.txt &&
			grep "## AverageGaps" temp_out_sgt_$filename.txt | awk '{print $3}' > avg_gaps_$filename.txt &&
			grep "## PercentageGappyColumns50" temp_out_sgt_$filename.txt | awk '{print $3}' > gappy_columns_50_$filename.txt &&
			grep "## PercentageGappyColumns80" temp_out_sgt_$filename.txt | awk '{print $3}' > gappy_columns_80_$filename.txt &&
			rm temp_out_sgt_$filename.txt) &
            $trimal_local -in $file -sident | grep "## AverageIdentity" | awk '{print $3}' > avg_seq_identity_$filename.txt &
            if [ -s $tree_file ]; then
                ref_tree=$(find $folder_path -name *.reference.nwk)
                (ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename.txt &&
                awk 'FNR == 3 {print $9}' temp_RF_distance_$filename.txt || 
                echo -1) > RF_distance_$filename.txt &
            fi
        else
            python3 $scripts/get_pre_alignment_statistics.py -i $file > temporal_max_min_col_$filename.txt
            cat temporal_max_min_col_$filename.txt | grep '## Minimum columns' | awk '{print $4}' > min_columns_$filename.txt
            cat temporal_max_min_col_$filename.txt | grep '## Maximum columns' | awk '{print $4}' > max_columns_$filename.txt
            rm temporal_max_min_col_$filename.txt
        fi 
    else
        echo "$file already processed" >> log.txt
    fi
    wait
    if [ -s $file.treefile ]; then
        rm temp_RF_distance_$filename.txt
    fi
    echo $filename >> filenames_to_write.txt
    end_file_time=$(date +%s.%N)
    file_time=$(echo "$end_file_time - $start_file_time" | bc)
    echo "$filename processed in $file_time seconds"
    echo "$filename processed in $file_time seconds" >> log.txt 
}


write_results() {
    file="filenames_to_write.txt"
    while read -r filename; do
		echo "writing $filename" >> log.txt &
		residue_taxon_problem=$(echo $filename | awk -F '.' '{print $1}')
		if [[ $filename == *.fa_err ]]; then
			echo "True" >> error_problem.txt &
			file_extension=".fa_err"
			cat number_sequences_"$residue_taxon_problem"_err.txt >> number_sequences.txt
		else
			echo "False" >> error_problem.txt &
			cat number_sequences_$residue_taxon_problem.txt >> number_sequences.txt &
		fi
		echo $filename | awk -F '_' '{print $1}' >> residue_type.txt &
		echo $filename | awk -F '_' '{print $2}' >> taxon.txt &
		echo $filename | awk -F '.' '{print $1}' | grep -o '[0-9]\+' >> problem_num.txt &
		if [[ $filename != *".seqs."* && $filename != *.fa_err  ]]; then
			MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
			MSA_filter_tool=""
			if [[ $MSA_tool == Guidance* ]]; then
				MSA_tool=$(echo $MSA_tool | awk -F 'Guidance' '{print $2}')
				echo $MSA_tool  >> MSA_tools.txt &
				MSA_filter_tool='Guidance'
				echo $MSA_filter_tool >> MSA_filter_tools.txt &
			else
				echo $MSA_tool >> MSA_tools.txt &
				MSA_filter_tool=$(echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}')
				echo $MSA_filter_tool >> MSA_filter_tools.txt &
			fi
			avg_gaps_diff=-1
			avg_gaps_diff_weighted=-1
			avg_seq_identity_diff=-1
			avg_seq_identity_diff_weighted=-1
			removed_columns=-1
			percent_conserved_columns=-1
			blocks_diff=-1
            RF_distance_diff=-1
			if [[ $MSA_filter_tool != "None" ]]; then
				if [[ $MSA_filter_tool == "Guidance" ]]; then
					if [[ $MSA_tool == "ClustalW" ]]; then
						avg_gaps_MSA_tool_file="avg_gaps_$residue_taxon_problem."$MSA_tool"2.fa.txt"
						avg_seq_identity_MSA_tool_file="avg_seq_identity_$residue_taxon_problem."$MSA_tool"2.fa.txt"
						number_columns_MSA_tool_file="number_columns_$residue_taxon_problem."$MSA_tool"2.fa.txt"
						number_blocks_MSA_tool_file="number_blocks_$residue_taxon_problem."$MSA_tool"2.fa.txt"
                        RF_distance_MSA_tool_file="RF_distance_$residue_taxon_problem."$MSA_tool"2.fa.txt"
					else
						avg_gaps_MSA_tool_file="avg_gaps_$residue_taxon_problem.$MSA_tool.fa.txt"
						avg_seq_identity_tool_file="avg_seq_identity_$residue_taxon_problem.$MSA_tool.fa.txt"
						number_columns_MSA_tool_file="number_columns_$residue_taxon_problem.$MSA_tool.fa.txt"
						number_blocks_MSA_tool_file="number_blocks_$residue_taxon_problem.$MSA_tool.fa.txt"
                        RF_distance_MSA_tool_file="RF_distance_$residue_taxon_problem.$MSA_tool.fa.txt"
					fi
				else
					file_suffix=$(echo $filename | awk -F '.' '{print $(NF-3)"."$(NF-2)"."$NF}')
					avg_gaps_MSA_tool_file="avg_gaps_$file_suffix.txt"
					avg_seq_identity_MSA_tool_file="avg_seq_identity_$file_suffix.txt"
					number_columns_MSA_tool_file="number_columns_$file_suffix.txt"
					number_blocks_MSA_tool_file="number_blocks_$file_suffix.txt"
                    RF_distance_MSA_tool_file="RF_distance_$file_suffix.txt"
				fi
				avg_gaps_MSA_filter_tool=$(cat avg_gaps_$filename.txt)
				avg_gaps_MSA_tool=$(cat $avg_gaps_MSA_tool_file)
				avg_gaps_diff=$(echo "$avg_gaps_MSA_tool - $avg_gaps_MSA_filter_tool" | bc)
				avg_seq_identity_MSA_filter_tool=$(cat avg_seq_identity_$filename.txt)
				avg_seq_identity_MSA_tool=$(cat $avg_seq_identity_MSA_tool_file)
				avg_seq_identity_diff=$(echo "$avg_seq_identity_MSA_tool - $avg_seq_identity_MSA_filter_tool" | bc)
				number_columns_filter_tool=$(cat number_columns_$filename.txt)
				number_columns_tool=$(cat $number_columns_MSA_tool_file)
				removed_columns=$(echo "$number_columns_tool - $number_columns_filter_tool" | bc)
				percent_conserved_columns=$(echo "$number_columns_filter_tool / $number_columns_tool" | bc -l)
				number_blocks_filter_tool=$(cat number_blocks_$filename.txt)
				number_blocks_tool=$(cat $number_blocks_MSA_tool_file)
				blocks_diff=$(echo "$number_blocks_tool - $number_blocks_filter_tool" | bc)
				RF_distance_filter_tool=$(cat RF_distance_$filename.txt)
				RF_distance_tool=$(cat $RF_distance_MSA_tool_file)
				RF_distance_diff=$(echo "$RF_distance_tool - $RF_distance_filter_tool" | bc)
				avg_gaps_diff_weighted=$(echo "$avg_gaps_diff * $percent_conserved_columns" | bc -l)
				avg_seq_identity_diff_weighted=$(echo "$avg_seq_identity_diff * $percent_conserved_columns" | bc -l)
			else
				number_columns_tool=$(cat number_columns_$filename.txt || echo -1)
			fi
			(echo $number_columns_tool >> number_columns_MSA.txt) &
			(echo $avg_gaps_diff >> avg_gaps_diff.txt) &
			(echo $avg_seq_identity_diff >> avg_seq_identity_diff.txt) &
			(echo $avg_gaps_diff_weighted >> avg_gaps_diff_weighted.txt) &
			(echo $avg_seq_identity_diff_weighted >> avg_seq_identity_diff_weighted.txt) &
			(echo $removed_columns >> removed_columns.txt) &
			(echo $percent_conserved_columns >> percent_conserved_columns.txt) &
			(echo $blocks_diff >> blocks_diff.txt) &
			if [[ -z $RF_distance_diff ]]; then
				(echo -100 >> RF_distance_diff.txt) &
			else
				(echo $RF_distance_diff >> RF_distance_diff.txt) &
			fi
			((cat number_blocks_$filename.txt || echo -1) >> number_blocks.txt) &
			((cat left_block_column_$filename.txt || echo -1) >> left_block_column.txt && rm left_block_column_$filename.txt) &
			((cat right_block_column_$filename.txt || echo -1) >> right_block_column.txt && rm right_block_column_$filename.txt) &
			((cat number_columns_$filename.txt || echo -1) >> min_columns.txt; (cat number_columns_$filename.txt || echo -1) >> max_columns.txt;
			(cat number_columns_$filename.txt || echo -1) >> number_columns.txt) &
			((cat avg_gaps_$filename.txt || echo -1) >> avg_gaps.txt) &
			((cat avg_seq_identity_$filename.txt || echo -1) >> avg_seq_identity.txt) &
			((cat gappy_columns_50_$filename.txt || echo -1) >> gappy_columns_50.txt && rm gappy_columns_50_$filename.txt) &
			((cat gappy_columns_80_$filename.txt || echo -1) >> gappy_columns_80.txt && rm gappy_columns_80_$filename.txt) &
            (cat RF_distance_$filename.txt || echo -1) >> RF_distance.txt
			rm temp_out_$filename.txt &
		else
			echo "Original" >> MSA_tools.txt &
			echo "Original" >> MSA_filter_tools.txt &
			(cat min_columns_$filename.txt >> min_columns.txt && rm min_columns_$filename.txt) &
			(cat max_columns_$filename.txt >> max_columns.txt && rm max_columns_$filename.txt) &
			echo -1 >> number_blocks.txt &
			echo -1 >> left_block_column.txt &
			echo -1 >> right_block_column.txt &
			echo -1 >> blocks_diff.txt &
			echo -1 >> number_columns.txt &
			echo -1 >> number_columns_MSA.txt &
			echo -1 >> removed_columns.txt &
			echo -1 >> percent_conserved_columns.txt &
			echo -1 >> avg_gaps.txt &
			echo -1 >> avg_gaps_diff.txt &
			echo -1 >> avg_gaps_diff_weighted.txt &
			echo -1 >> avg_seq_identity.txt &
			echo -1 >> avg_seq_identity_diff.txt &
			echo -1 >> avg_seq_identity_diff_weighted.txt &
			echo -1 >> RF_distance.txt &
            echo -1 >> RF_distance_diff.txt &
			echo -1 >> gappy_columns_50.txt &
			echo -1 >> gappy_columns_80.txt &
		fi
		wait
        echo "Finished with $filename"
        echo "Finished with $filename" >> log.txt
    done < $file
	# skip if no file written
	if [ -s filenames_to_write.txt ]; then
		rm avg_gaps_*problem*.txt
		rm avg_seq_identity_*problem*.txt
		rm number_columns_*problem*.txt
		rm number_blocks_*.txt
		rm RF_distance_*problem*.txt
    fi

    > filenames_to_write.txt
}


counter=0
line_number=1
for residue_type in $dataset/$res_type*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/Fun*
        do
	    	residue_taxon_problem=""
            for problem in $taxon/problem093*
            do
                problem_name=$(basename $problem | awk -F '_' '{print $1}')
				residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
				#residue_taxon_problem_to_ignore=$(sed "${line_number}q;d" $problems_to_ignore)

				if grep -Fxq $residue_taxon_problem $problems_to_ignore; then
					((line_number++))
					echo "Ignored $residue_taxon_problem" >> log.txt
				else
					if [[ $problem == *"problem"*"_err" ]]; then
						num_seq=$(grep -c ">" "$problem/$problem_name.seqs.fa_err")
						echo $num_seq > number_sequences_$residue_taxon_problem.txt
					else
						num_seq=$(grep -c ">" "$problem/$problem_name.seqs.fa")
						echo $num_seq > number_sequences_$residue_taxon_problem.txt
						for file in $problem/*.fa
						do
							task "$file"  &
						done
					fi
					((counter++))
				fi
				
				if [ $counter -eq 2 ]; then
					wait
					counter=0
					write_results
					rm number_sequences_*.txt
				fi
            done
		done
    fi
done

wait
write_results
rm number_sequences_*.txt

end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
