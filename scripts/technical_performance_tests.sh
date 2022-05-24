#!/bin/bash

file=$1

start_time=$(date +%s.%N)
trimal_local="trimal/bin/trimal"
time_command="/usr/bin/time"

methods=("None" "strict" "strictplus" "gappyout")

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
        (cat ${exec_times_filename}.txt >> ${exec_times_filename}.txt && rm ${exec_times_filename}$filename.txt)
        (cat ${exec_times_filename}.txt >> ${exec_times_filename}.txt && rm ${exec_times_filename}_$filename.txt)
        (cat alignment_sizes_$filename.txt >> alignment_sizes.txt && rm alignment_sizes_$filename.txt)
        (cat alignment_seqs_$filename.txt >> alignment_seqs.txt && rm alignment_seqs_$filename.txt)
        (cat alignment_cols_$filename.txt >> alignment_cols.txt && rm alignment_cols_$filename.txt)
        (cat methods_$filename.txt >> methods.txt && rm methods_$filename.txt)
        wait
    done < $file
    > filenames_to_write.txt
}

test_performance() {
    file=$1
    filename=$(echo $file | awk -F '/' '{print $NF}')
    alignment_size=$(stat --printf="%s\n" $file)
    alignment_seqs=$(echo $filename | awk -F '.' '{print $(NF-2)}')
    alignment_cols=$(echo $filename | awk -F '.' '{print $(NF-1)}')

    for method in ${methods[@]}
    do
        if [[ $method == "None" ]]; then
            exec_times_filename="exec_times_$filename.txt"
            { $time_command -v $trimal_local -in $file >/dev/null ; } 2> $exec_times_filename
        else
            exec_times_filename="exec_times_${method}_$filename.txt"
            { $time_command -v $trimal_local -in $file -$method >/dev/null ; } 2> $exec_times_filename
        fi

        user_time=$(grep "User time" $exec_times_filename  | awk -F ':' '{print $NF}')
        system_time=$(grep "System time" $exec_times_filename  | awk -F ':' '{print $NF}')
        percent_cpu=$(grep "Percent of CPU" $exec_times_filename  | awk -F ':' '{print $NF}' | awk -F '%' '{print $(NF-1)}')
        max_resident_set_size=$(grep "Maximum resident" $exec_times_filename  | awk -F ':' '{print $NF}')
        exit_status=$(grep "Exit status" $exec_times_filename | awk -F ':' '{print $NF}')

        echo "$file,$method,$alignment_size,$alignment_seqs,$alignment_cols,$user_time,$system_time,$percent_cpu,$max_resident_set_size,$exit_status"
    done
}

# TODO: Repeat 100 times
test_performance $file
