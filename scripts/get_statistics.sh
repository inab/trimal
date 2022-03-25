#!/bin/bash

start_time=$(date +%s.%N)
scripts="../trimal/scripts"
dataset="../dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
test_working_files="test_working_files"
trimal_local="../trimal/bin/trimal"
if [ ! -d "$test_working_files" ]; then
    mkdir $test_working_files
fi
cd $test_working_files

> blocks_outputs.txt
> left_block_column.txt
> right_block_column.txt
> number_sequences.txt
> number_columns.txt
> max_columns.txt
> min_columns.txt
> avg_gaps.txt
> avg_gaps_diff.txt
> avg_seq_identity.txt
> gappy_columns_50.txt
> gappy_columns_80.txt
> RF_distance.txt
> residue_type.txt
> taxon.txt
> problem_num.txt
> error_problem.txt
> MSA_tools.txt
> MSA_filter_tools.txt
> log.txt
> log_processing.txt
> log_writing.txt
> filenames_to_write.txt


task(){
    file=$1
    start_file_time=$(date +%s.%N)
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    if [ ! -s temp_out_$filename.txt ]; then
        folder_path=$(dirname $file)
        if [[ $file == *".fa" ]]; then
            filename_no_ext=$(basename $file .fa)
            tree_file="$folder_path/$filename_no_ext.nwk"
        fi
        if [[ $filename != *".seqs."* ]]; then
            $trimal_local -in $file -sgc > temp_out_$filename.txt
            (python3 $scripts/set_manual_boundaries.py -i temp_out_$filename.txt --inner_blocks --total_blocks --min_gapscore_allowed 1 > temp_blocks_outputs_$filename.txt && 
			grep '## Blocks' temp_blocks_outputs_$filename.txt | awk '{print $3}' > blocks_outputs_$filename.txt &&
			grep '## Left column' temp_blocks_outputs_$filename.txt | awk '{print $4}' > left_block_column_$filename.txt &&
			grep '## Right column' temp_blocks_outputs_$filename.txt | awk '{print $4}' > right_block_column_$filename.txt &&
			rm temp_blocks_outputs_$filename.txt) &
            tail -n 1 temp_out_$filename.txt | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' > number_columns_$filename.txt &
			($trimal_local -in $file -sgt > temp_out_sgt_$filename.txt &&
			grep "## AverageGaps" temp_out_sgt_$filename.txt | awk '{print $3}' > avg_gaps_$filename.txt &&
			grep "## PercentageGappyColumns50" temp_out_sgt_$filename.txt | awk '{print $3}' > gappy_columns_50_$filename.txt &&
			grep "## PercentageGappyColumns80" temp_out_sgt_$filename.txt | awk '{print $3}' > gappy_columns_80_$filename.txt &&
			rm temp_out_sgt_$filename.txt) &
            $trimal_local -in $file -sident | grep "## AverageIdentity" | awk '{print $3}' > avg_seq_identity_$filename.txt &
            if [ -s $tree_file ]; then
                ref_tree=$(find $folder_path -name *.reference.nwk)
                ete3 compare -t $tree_file -r $ref_tree | awk 'FNR == 3 {print $9}' > RF_distance_$filename.txt &
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
    echo $filename >> filenames_to_write.txt
    end_file_time=$(date +%s.%N)
    file_time=$(echo "$end_file_time - $start_file_time" | bc)
    echo "$filename processed in $file_time seconds" >> log_processing.txt 
}


write_results() {
    file="filenames_to_write.txt"
    while read -r filename; do
		echo "writing $filename" >> log.txt &
		residue_taxon_problem=$(echo $filename | awk -F '.' '{print $1}')
		if [[ $filename == *.fa_err ]]; then
			echo "True" >> error_problem.txt &
			cat number_sequences_"$residue_taxon_problem"_err.txt >> number_sequences.txt &
		else
			echo "False" >> error_problem.txt &
			cat number_sequences_$residue_taxon_problem.txt >> number_sequences.txt &
		fi
		echo $filename | awk -F '_' '{print $1}' >> residue_type.txt &
		echo $filename | awk -F '_' '{print $2}' >> taxon.txt &
		echo $filename | awk -F '.' '{print $1}' | grep -o '[0-9]\+' >> problem_num.txt &
		if [[ $filename != *".seqs."* ]]; then
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
			if [[ $MSA_filter_tool != "None" ]]; then
				if [[ $MSA_filter_tool == "Guidance" ]]; then
					if [[ $MSA_tool == "ClustalW" ]]; then
						avg_gaps_MSA_tool_file="avg_gaps_$residue_taxon_problem."$MSA_tool"2.fa.txt"
					else
						avg_gaps_MSA_tool_file="avg_gaps_$residue_taxon_problem.$MSA_tool.fa.txt"
					fi

					echo $filename
					#echo "avg_gaps_MSA_tool=$avg_gaps_MSA_tool"
					avg_gaps_MSA_tool=$(cat "$avg_gaps_MSA_tool_file")
					avg_gaps_MSA_filter_tool=$(cat avg_gaps_$filename.txt)

				else
					avg_gaps_MSA_filter_tool=$(cat avg_gaps_$filename.txt)
					avg_gaps_MSA_tool_file=$(echo $filename | awk -F '.' '{print $(NF-3)"."$(NF-2)"."$NF}')
					avg_gaps_MSA_tool_file="avg_gaps_$avg_gaps_MSA_tool_file.txt"
					echo $filename
					#echo "avg_gaps_MSA_tool=$avg_gaps_MSA_tool"
					avg_gaps_MSA_tool=$(cat "$avg_gaps_MSA_tool_file")
					avg_gaps_MSA_filter_tool=$(cat avg_gaps_$filename.txt)
				fi
				avg_gaps_diff=$(echo "$avg_gaps_MSA_tool - $avg_gaps_MSA_filter_tool" | bc)
			fi
			(echo $avg_gaps_diff >> avg_gaps_diff.txt) &
			(cat blocks_outputs_$filename.txt >> blocks_outputs.txt && rm blocks_outputs_$filename.txt) &
			(cat left_block_column_$filename.txt >> left_block_column.txt && rm left_block_column_$filename.txt) &
			(cat right_block_column_$filename.txt >> right_block_column.txt && rm right_block_column_$filename.txt) &
			(cat number_columns_$filename.txt >> min_columns.txt; cat number_columns_$filename.txt >> max_columns.txt;
			cat number_columns_$filename.txt >> number_columns.txt && rm number_columns_$filename.txt) &
			(cat avg_gaps_$filename.txt >> avg_gaps.txt) &
			(cat avg_seq_identity_$filename.txt >> avg_seq_identity.txt && rm avg_seq_identity_$filename.txt) &
			(cat gappy_columns_50_$filename.txt >> gappy_columns_50.txt && rm gappy_columns_50_$filename.txt) &
			(cat gappy_columns_80_$filename.txt >> gappy_columns_80.txt && rm gappy_columns_80_$filename.txt) &
			if [ -s RF_distance_$filename.txt ]; then
				cat RF_distance_$filename.txt >> RF_distance.txt && rm RF_distance_$filename.txt
			else
				echo -1 >> RF_distance.txt &
			fi
			rm temp_out_$filename.txt &
		else
			echo "Original" >> MSA_tools.txt &
			echo "Original" >> MSA_filter_tools.txt &
			(cat min_columns_$filename.txt >> min_columns.txt && rm min_columns_$filename.txt) &
			(cat max_columns_$filename.txt >> max_columns.txt && rm max_columns_$filename.txt) &
			echo -1 >> blocks_outputs.txt &
			echo -1 >> left_block_column.txt &
			echo -1 >> right_block_column.txt &
			echo -1 >> number_columns.txt &
			echo -1 >> avg_gaps.txt &
			echo -1 >> avg_gaps_diff.txt &
			echo -1 >> avg_seq_identity.txt &
			echo -1 >> RF_distance.txt &
			echo -1 >> gappy_columns_50.txt &
			echo -1 >> gappy_columns_80.txt &
		fi
		echo "Finished with $filename" >> log_writing.txt &
		wait
    done < $file
	rm avg_gaps_*problem*.txt
    > filenames_to_write.txt
}


counter=0
for residue_type in $dataset/*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/Bac*
        do
	    	residue_taxon_problem=""
            for problem in $taxon/problem0001*
            do
                problem_name=$(basename $problem | awk -F '_' '{print $1}')
				residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
                if [[ $problem == *"problem"*"_err" ]]; then
                    num_seq=$(grep ">" "$problem/$problem_name.seqs.fa_err" | wc -l)
                    echo $num_seq > number_sequences_$residue_taxon_problem.txt
		    		for file in $problem/*.fa_err
                    do
                        task "$file"  &
                    done
                else
                    num_seq=$(grep ">" "$problem/$problem_name.seqs.fa" | wc -l)
		    		echo $num_seq > number_sequences_$residue_taxon_problem.txt
                    for file in $problem/*.fa
                    do
                        task "$file"  &
                    done
                fi

				((counter++))
				if [ $counter -eq 100 ]; then
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
