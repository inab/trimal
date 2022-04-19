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
> methods.txt


methods=("None" "strict" "strictplus" "gappyout" "automated1")

parse_results() {
    method=$1
    if [[ $method == "None" ]]; then
        exec_times_filename="exec_times.txt"
    else
        exec_times_filename="exec_times_$method.txt"
    fi
    grep "Command being" $exec_times_filename  | awk -F '/' '{print $NF}' | awk -F '"' '{print $(NF-1)}' > alignment_filenames.txt
    grep "User time" $exec_times_filename  | awk -F ':' '{print $NF}' > user_times.txt
    grep "System time" $exec_times_filename  | awk -F ':' '{print $NF}' > system_times.txt
    grep "Percent of CPU" $exec_times_filename  | awk -F ':' '{print $NF}' | awk -F '%' '{print $(NF-1)}' > percent_cpu.txt
    grep "Maximum resident" $exec_times_filename  | awk -F ':' '{print $NF}' > max_resident_set_size.txt
    grep "Exit status" $exec_times_filename | awk -F ':' '{print $NF}' > exit_status.txt
}

write_results() {
    file="filenames_to_write.txt"
    method=$1
    if [[ $method == "None" ]]; then
        exec_times_filename="exec_times_$filename.txt"
    else
        exec_times_filename="exec_times_${method}_$filename.txt"
    fi
    while read -r filename; do
        (cat ${exec_times_filename}.txt >> ${exec_times_filename}.txt && rm ${exec_times_filename}$filename.txt) &
        (cat ${exec_times_filename}.txt >> ${exec_times_filename}.txt && rm ${exec_times_filename}_$filename.txt) &
        (cat alignment_sizes_$filename.txt >> alignment_sizes.txt && rm alignment_sizes_$filename.txt) &
        (cat alignment_seqs_$filename.txt >> alignment_seqs.txt && rm alignment_seqs_$filename.txt) &
        (cat alignment_cols_$filename.txt >> alignment_cols.txt && rm alignment_cols_$filename.txt) &
        (cat methods_$filename.txt >> methods.txt && rm methods_$filename.txt) &
        wait
    done < $file
    > filenames_to_write.txt
}

test_performance() {
    file=$1
    filename=$(echo $file | awk -F '/' '{print $NF}')
    alignment_size=$(stat --printf="%s\n" $file) &
    alignment_seqs=$(grep -c ">" $file) &
    alignment_cols=$($trimal_local -in $file -sgc | tail -n 1 | grep -E "^\s*[[:digit:]]{1,}" -o | awk '{print $1 + 1}' || echo -1) &
    echo $method >> methods_$filename.txt
    for method in ${methods[@]}
    do
        if [[ $method == "None" ]]; then
            exec_times_filename="exec_times_$filename.txt"
            { $time_command -v $trimal_local -in $file >/dev/null ; } 2>> $exec_times_filename.txt &
        else
            exec_times_filename="exec_times_${method}_$filename.txt"
            { $time_command -v $trimal_local -in $file -$method >/dev/null ; } 2>> $exec_times_filename.txt &
        fi
    done

    wait
    for method in ${methods[@]}
    do
        echo $alignment_size >> alignment_sizes_$filename.txt
        echo $alignment_seqs >> alignment_seqs_$filename.txt
        echo $alignment_cols >> alignment_cols_$filename.txt
        echo $method >> methods.txt
    done

    echo $filename >> filenames_to_write.txt
}

counter=0
for file in $dataset/*
do
    test_performance $file &
    ((counter++))
    if [[ $counter -eq 3 ]]; then
        wait
        counter=0
        write_results
        break
    fi
done
wait
write_results $method

for method in ${methods[@]}
do
    parse_results $method
done


end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
