#!/bin/bash

start_time=$(date +%s.%N)
scripts="../trimal/scripts"
dataset="/home/shared/dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
test_working_files="test_working_files"
trimal_local="../trimal/source/trimal"
if [ ! -d "$test_working_files" ]; then
    mkdir $test_working_files
fi
cd $test_working_files

> blocks_outputs.txt
> number_sequences.txt
> number_columns.txt
> max_columns.txt
> min_columns.txt
> avg_gaps.txt
> avg_seq_identity.txt
> RF_distance.txt
> residue_type.txt
> taxon.txt
> problem_num.txt
> error_problem.txt
> MSA_tools.txt
> MSA_filter_tools.txt
> log_processing.txt
> log_writing.txt

task(){
    file=$1
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    if [ ! -s temporal_out_$filename.txt ]; then
        folder_path=$(dirname $file)
        if [[ $file == *".fa" ]]; then
            filename_no_ext=$(basename $file .fa)
            tree_file="$folder_path/$filename_no_ext.nwk"
        fi
        if [[ $filename != *".seqs."* ]]; then
            $trimal_local -in $file -sgc > temporal_out_$filename.txt
            python3 $scripts/set_manual_boundaries.py -i temporal_out_$filename.txt --inner_blocks --total_blocks --min_gapscore_allowed 1 > blocks_outputs_$filename.txt &
            cat temporal_out_$filename.txt | tail -n 1 | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' > number_columns_$filename.txt &
            #{
                #trimal -in $file -sst -sident -soverlap > temporal_sst_sident_soverlap_$filename.txt;
                #grep '^[[:digit:]]*\s' temporal_sst_sident_soverlap_$filename.txt | awk '{sum+=$1;} END{print sum}' > number_columns_$filename.txt;
                #grep '^[[:digit:]]*\s' temporal_sst_sident_soverlap_$filename.txt | awk 'NR==1' | awk '{if ($5 == 1) {print $1} else {print 0}}' > number_identical_columns_$filename.txt;
                #grep "## AverageIdentity" temporal_sst_sident_soverlap_$filename.txt | awk '{print $3}' > avg_seq_identity_$filename.txt;
                #grep "## AverageOverlap" temporal_sst_sident_soverlap_$filename.txt | awk '{print $3}' > avg_seq_overlap_$filename.txt;
                #rm temporal_sst_sident_soverlap_$filename.txt;
            #}
            $trimal_local -in $file -sgt | grep "## AverageGaps" | awk '{print $3}' > avg_gaps_$filename.txt &
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
        echo "$file already processed"
    fi
    wait
}


N=5
for residue_type in $dataset/*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem*
            do
                start_problem_time=$(date +%s.%N)
                problem_name=$(basename $problem | awk -F '_' '{print $1}')
                problem_num=$(echo $problem_name | grep -o '[0-9]\+')
                if [[ $problem == *"problem"*"_err" ]]; then
                    num_seq=$(grep ">" "$problem/$problem_name.seqs.fa_err" | wc -l)
                    for file in $problem/*.fa_err
                    do
                        ((i=i%N)); ((i++==0)) && wait
                        echo $num_seq >> number_sequences.txt
                        task "$file"  &
                    done
                else
                    num_seq=$(grep ">" "$problem/$problem_name.seqs.fa" | wc -l)
                    for file in $problem/*.fa
                    do
                        ((i=i%N)); ((i++==0)) && wait
                        echo $num_seq >> number_sequences.txt
                        task "$file"  &
                    done
                fi

                end_problem_time=$(date +%s.%N)
                problem_time=$(echo "$end_problem_time - $start_problem_time" | bc)
                echo "$problem processed in $problem_time seconds" >> log_processing.txt
            done
        done
    fi
done
wait


for residue_type in $dataset/*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem*
            do
                file_format=".fa"
                if [[ $problem == *"problem"*"_err" ]]; then
                    file_format=".fa_err"
                fi
                for file in $problem/*$file_format
                do
                    if [[ $problem == *"problem"*"_err" ]]; then
                        echo "True" >> error_problem.txt
                    else
                        echo "False" >> error_problem.txt
                    fi
                    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
                    echo $residue_type | awk -F '/' '{print $(NF)}' >> residue_type.txt
                    echo $taxon | awk -F '/' '{print $(NF)}' >> taxon.txt 
                    basename $problem | grep -o '[0-9]\+' >> problem_num.txt
                    if [[ $filename != *".seqs."* ]]; then
                        MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
                        if  [[ $MSA_tool == Guidance* ]]; then
                            echo $MSA_tool | awk -F 'Guidance' '{print $2}' >> MSA_tools.txt
                            echo 'Guidance' >> MSA_filter_tools.txt
                        else
                            echo $MSA_tool >> MSA_tools.txt
                            echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}' >> MSA_filter_tools.txt
                        fi
                        cat blocks_outputs_$filename.txt >> blocks_outputs.txt && rm blocks_outputs_$filename.txt
                        cat number_columns_$filename.txt >> min_columns.txt
                        cat number_columns_$filename.txt >> max_columns.txt
                        cat number_columns_$filename.txt >> number_columns.txt && rm number_columns_$filename.txt
                        cat avg_gaps_$filename.txt >> avg_gaps.txt && rm avg_gaps_$filename.txt
                        cat avg_seq_identity_$filename.txt >> avg_seq_identity.txt && rm avg_seq_identity_$filename.txt
                        if [ -s RF_distance_$filename.txt ]; then
                            cat RF_distance_$filename.txt >> RF_distance.txt && rm RF_distance_$filename.txt
                        else
                            echo -1 >> RF_distance.txt
                        fi
                        rm temporal_out_$filename.txt
                    else
                        echo "Original" >> MSA_tools.txt
                        echo "Original" >> MSA_filter_tools.txt
                        cat min_columns_$filename.txt >> min_columns.txt && rm min_columns_$filename.txt
                        cat max_columns_$filename.txt >> max_columns.txt && rm max_columns_$filename.txt
                        echo -1 >> blocks_outputs.txt
                        echo -1 >> number_columns.txt
                        echo -1 >> avg_gaps.txt
                        echo -1 >> avg_seq_identity.txt
                        echo -1 >> RF_distance.txt
                    fi
                    echo "Finished with $filename" >> log_writing.txt
                done
            done
        done
    fi
done


end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
echo $diff