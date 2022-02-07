#!/bin/bash

start_time=$(date +%s.%N)
scripts="../trimal/scripts"
dataset="/home/shared/dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched/AA"
test_working_files="test_working_files"
trimal_local="../trimal/source/trimal"
if [ ! -d "$test_working_files" ]; then
    mkdir $test_working_files
fi
cd $test_working_files

task(){
    start_file_time=$(date +%s.%N)
    file=$1
    filename=$(echo $file | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
    $trimal_local -in $file -sgc > temporal_out_$filename.txt
    if [ -s temporal_out_$filename.txt ]; then
        python3 $scripts/set_manual_boundaries.py -i temporal_out_$filename.txt --inner_blocks --total_blocks --min_gapscore_allowed 1 > blocks_outputs_$filename.txt &
        grep ">" $file | wc -l > number_sequences_$filename.txt &
        $trimal_local -in $file -sgc | tail -n 1 | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' > number_columns_$filename.txt &
        #{
            #trimal -in $file -sst -sident -soverlap > temporal_sst_sident_soverlap_$filename.txt;
            #grep '^[[:digit:]]*\s' temporal_sst_sident_soverlap_$filename.txt | awk '{sum+=$1;} END{print sum}' > number_columns_$filename.txt;
            #grep '^[[:digit:]]*\s' temporal_sst_sident_soverlap_$filename.txt | awk 'NR==1' | awk '{if ($5 == 1) {print $1} else {print 0}}' > number_identical_columns_$filename.txt;
            #grep "## AverageIdentity" temporal_sst_sident_soverlap_$filename.txt | awk '{print $3}' > avg_seq_identity_$filename.txt;
            #grep "## AverageOverlap" temporal_sst_sident_soverlap_$filename.txt | awk '{print $3}' > avg_seq_overlap_$filename.txt;
            #rm temporal_sst_sident_soverlap_$filename.txt;
        #}
        $trimal_local -in $file -sgt | grep "## AverageGaps" | awk '{print $3}' > avg_gaps_$filename.txt &
        #$trimal_local -in $file -sident | grep "## AverageIdentity" | awk '{print $3}' > avg_seq_identity_$filename.txt &
    fi
    wait
    end_file_time=$(date +%s.%N)
    file_time=$(echo "$end_file_time - $start_file_time" | bc)
    filename="${file##*/}"
    #echo "$filename took $file_time  seconds"
}

for taxon in $dataset/*
do
    for problem in $taxon/problem
    do
        for file in $problem/*.fa
        do
            task "$file" &
        done
    done
done
wait


> blocks_outputs.txt
> number_sequences.txt
> number_columns.txt
> avg_gaps.txt
> taxon.txt
> problem_num.txt
> MSA_tools.txt
> MSA_filter_tools.txt

for taxon in $dataset/*
do
    for problem in $taxon/problem*
    do
        for file in $problem/*.fa
        do
            filename=$(echo $file | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
            echo $filename
            if [ -s temporal_out_$filename.txt ]; then
                echo $taxon | awk -F '/' '{print $(NF)}' >> taxon.txt 
                echo $problem | awk -F '/' '{print $(NF)}' >> problem_num.txt
                MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
                if  [[ $MSA_tool == Guidance* ]]; then
                    MSA_tool=$(echo $MSA_tool | awk -F 'Guidance' '{print $2}')
                fi
                echo $MSA_tool >> MSA_tools.txt
                echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}' >> MSA_filter_tools.txt
                cat blocks_outputs_$filename.txt >> blocks_outputs.txt && rm blocks_outputs_$filename.txt
                cat number_sequences_$filename.txt >> number_sequences.txt && rm number_sequences_$filename.txt
                cat number_columns_$filename.txt >> number_columns.txt && rm number_columns_$filename.txt
                cat avg_gaps_$filename.txt >> avg_gaps.txt && rm avg_gaps_$filename.txt
            fi
            rm temporal_out_$filename.txt
        done
    done
done


end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff