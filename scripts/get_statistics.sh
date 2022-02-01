#!/bin/bash

start_time=$(date +%s.%N)
scripts="/trimal/scripts/"
dataset="/trimal/dataset/"
test_working_files="test_working_files"
cd $test_working_files

task(){
    start_file_time=$(date +%s.%N)
    file=$1
    #echo $file
    filename="${file##*/}"
    trimal -in $file -sgc > temporal_out_$filename.txt
    if [ -s temporal_out_$filename.txt ]; then
        python3 $scripts/set_manual_boundaries.py -i temporal_out_$filename.txt --inner_blocks --total_blocks --min_gapscore_allowed 1 > blocks_outputs_$filename.txt &
        grep ">" $file | wc -l > number_sequences_$filename.txt &
        trimal -in $file -sgc | tail -n 1 | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' > number_columns_$filename.txt &
        trimal -in $file -sst | grep '[[:digit:]]' | awk 'NR==1' | awk '{if ($5 == 1) {print $1} else {print 0}}' > number_identical_columns_$filename.txt &
        #{
            #trimal -in $file -sst -sident -soverlap > temporal_sst_sident_soverlap_$filename.txt;
            #grep '^[[:digit:]]*\s' temporal_sst_sident_soverlap_$filename.txt | awk '{sum+=$1;} END{print sum}' > number_columns_$filename.txt;
            #grep '^[[:digit:]]*\s' temporal_sst_sident_soverlap_$filename.txt | awk 'NR==1' | awk '{if ($5 == 1) {print $1} else {print 0}}' > number_identical_columns_$filename.txt;
            #grep "## AverageIdentity" temporal_sst_sident_soverlap_$filename.txt | awk '{print $3}' > avg_seq_identity_$filename.txt;
            #grep "## AverageOverlap" temporal_sst_sident_soverlap_$filename.txt | awk '{print $3}' > avg_seq_overlap_$filename.txt;
            #rm temporal_sst_sident_soverlap_$filename.txt;
        #}
        trimal -in $file -sident | grep "## AverageIdentity" | awk '{print $3}' > avg_seq_identity_$filename.txt &
    fi
    wait
    end_file_time=$(date +%s.%N)
    file_time=$(echo "$end_file_time - $start_file_time" | bc)
    echo "$file took $file_time  seconds"
}

for file in $dataset*.fasta;
do
    task "$file" &
done
wait


> blocks_outputs.txt
> number_sequences.txt
> number_columns.txt
> avg_seq_identity.txt
> number_identical_columns.txt

for file in $dataset*.fasta;
do
    filename="${file##*/}"
    cat blocks_outputs_$filename.txt >> blocks_outputs.txt && rm blocks_outputs_$filename.txt
    cat number_sequences_$filename.txt >> number_sequences.txt && rm number_sequences_$filename.txt
    cat number_columns_$filename.txt >> number_columns.txt && rm number_columns_$filename.txt
    cat avg_seq_identity_$filename.txt >> avg_seq_identity.txt && rm avg_seq_identity_$filename.txt
    cat number_identical_columns_$filename.txt >> number_identical_columns.txt && rm number_identical_columns_$filename.txt
    rm temporal_out_$filename.txt
done


end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff