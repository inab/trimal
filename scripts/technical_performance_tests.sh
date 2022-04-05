#!/bin/bash

start_time=$(date +%s.%N)
dataset="../cedric_datasets/homfam/alignments"
technical_performance_results="technical_performance_results"
trimal_local="../trimal/bin/trimal"
time_command="/usr/bin/time"
if [ ! -d "$technical_performance_results" ]; then
    mkdir $technical_performance_results
fi
cd $technical_performance_results

> filenames_to_write.txt
> exec_times.txt
> alignment_sizes.txt
> alignment_seqs.txt
> alignment_cols.txt


parse_results() {
    grep "Command being" exec_times.txt | awk -F '/' '{print $NF}' | awk -F '"' '{print $(NF-1)}' > alignment_filenames.txt
    grep "User time" exec_times.txt | awk -F ':' '{print $NF}' > user_times.txt
    grep "System time" exec_times.txt | awk -F ':' '{print $NF}' > system_times.txt
    grep "Percent of CPU" exec_times.txt | awk -F ':' '{print $NF}' | awk -F '%' '{print $(NF-1)}' > percent_cpu.txt
    grep "Maximum resident" exec_times.txt | awk -F ':' '{print $NF}' > max_resident_set_size.txt
    grep "Exit status" exec_times.txt | awk -F ':' '{print $NF}' > exit_status.txt
}

write_results() {
    file="filenames_to_write.txt"
    while read -r filename; do
        (cat exec_times_$filename.txt >> exec_times.txt && rm exec_times_$filename.txt) &
        (cat alignment_sizes_$filename.txt >> alignment_sizes.txt && rm alignment_sizes_$filename.txt) &
        (cat alignment_seqs_$filename.txt >> alignment_seqs.txt && rm alignment_seqs_$filename.txt) &
        (cat alignment_cols_$filename.txt >> alignment_cols.txt && rm alignment_cols_$filename.txt) &
        wait
    done < $file
    > filenames_to_write.txt
}

test_performance() {
    file=$1
    filename=$(echo $file | awk -F '/' '{print $NF}')
    { $time_command -v $trimal_local -in $file >/dev/null ; } 2>> exec_times_$filename.txt &
    stat --printf="%s\n" $file >> alignment_sizes_$filename.txt &
    grep -c ">" $file >> alignment_seqs_$filename.txt &
    ($trimal_local -in $file -sgc | tail -n 1 | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' || echo -1) >> alignment_cols_$filename.txt
    wait
    echo $filename >> filenames_to_write.txt
}

counter=0
for file in $dataset/*
do
    test_performance $file &
    ((counter++))
    if [[ $counter -eq 48 ]]; then
        wait
        counter=0
        write_results
    fi
done
wait
write_results
parse_results


end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
