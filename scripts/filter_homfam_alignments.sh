#!/bin/bash

file=$1

dataset_ref_alignments="cedric_datasets/refs"
trimal_local="trimal/bin/trimal"
filter_reference_sequences_script="trimal/scripts/filter_homfam_alignments.py"


compute_seqs_stats() {
    file=$1
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    python3 $filter_reference_sequences_script -i $file -r $ref_alignment_name.ref
    $trimal_local -in ${file}_filtered_refs -noallgaps > "${file}_no_all_gaps_cols"
    clean_num_seqs=$(grep -c ">" ${file}_no_all_gaps_cols)

    trimmed_num_seqs_gappyout=$($trimal_local -in ${file}_no_all_gaps_cols -gappyout | grep -c ">")
    trimmed_num_seqs_strictplus=$($trimal_local -in ${file}_no_all_gaps_cols -strictplus | grep -c ">")

    rm ${file}_filtered_refs
    rm ${file}_no_all_gaps_cols

    echo "$file,$clean_num_seqs,$trimmed_num_seqs_gappyout,$trimmed_num_seqs_strictplus" 
}

compute_sum_of_pairs_stats() {
    file=$1
    original_num_seqs=$(grep -c ">" $file)
    original_num_cols=$(while read line; do echo $line; [ -z $line ] && break; done <<<  "$(awk -F '>' '{print $1}' ${file} | tail -n +2)" | paste -sd '' | wc -c)
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    python3 $filter_reference_sequences_script -i $file -r $ref_alignment_name.ref
    $trimal_local -in ${file}_filtered_refs -noallgaps > "${file}_no_all_gaps_cols"
    clean_num_seqs=$(grep -c ">" ${file}_no_all_gaps_cols)
    #clean_num_cols=$(while read line; do echo $line; [ -z $line ] && break; done <<<  "$(awk -F '>' '{print $1}' ${file}_no_all_gaps_cols | tail -n +2)" | paste -sd '' | wc -c)

    filename=$(echo $file | awk -F '/' '{print $NF}')
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    if [ ! -s ref_alignment_SoP_$ref.txt ]; then
        echo $ref_alignment_name.ref > ref_alignment_SoP_$ref.txt
    fi
    
    $trimal_local -sfc -compareset ref_alignment_SoP_$ref.txt -forceselect ${file}_no_all_gaps_cols > SoP_output_${filename}.txt
    tail -n +9 SoP_output_${filename}.txt | awk '{print $2}' > SoP_by_column_${filename}.txt
    clean_num_cols=$(wc -l SoP_by_column_${filename}.txt | awk '{ print $1 }')
    if [[ $(echo "$clean_num_cols % 2" | bc) == '0' ]]; then
        median_pos_first=$(echo "$clean_num_cols / 2" | bc)
        median_pos_second=$(echo "$median_pos_first + 1" | bc)
        median_first=$(sed "${median_pos_first}q;d" SoP_by_column_${filename}.txt)
        median_second=$(sed "${median_pos_second}q;d" SoP_by_column_${filename}.txt)
        median_alignment=$(echo "($median_first + $median_second) / 2" | bc)
    else
        median_pos=$(echo "$clean_num_cols / 2 + 1" | bc)
        median_alignment=$(sed "${median_pos}q;d" SoP_by_column_${filename}.txt)
    fi

    min_SoP_alignment=$(sort SoP_by_column_${filename}.txt | head -n 1)
    max_SoP_alignment=$(sort SoP_by_column_${filename}.txt | tail -n 1)
    SoP_by_column_sum=$(paste -sd+ SoP_by_column_${filename}.txt | bc)
    alignment_SoP=$(echo "$SoP_by_column_sum / $clean_num_cols" | bc -l)
    std_alignment=$(awk -v avg_SoP=$alignment_SoP -v size=$clean_num_cols '{n += ($1 - avg_SoP) ** 2}; END{print sqrt(n/size)}' SoP_by_column_${filename}.txt)


    trimmed_column_positions_gappyout=$($trimal_local -in ${file}_no_all_gaps_cols -gappyout -colnumbering | grep "#ColumnsMap" | awk '{print substr($0,13)}')
    if [ -n "$trimmed_column_positions_gappyout" ]; then
        > SoP_by_column_${filename}_trimmed_gappyout.txt
        # -colnumbering returns columns starting by 0 and -sfc starting by 1
        for column_position in ${trimmed_column_positions_gappyout//,/};
        do
            col_pos=$((column_position + 1))
            sed "${col_pos}q;d" SoP_by_column_$filename.txt >> SoP_by_column_${filename}_trimmed_gappyout.txt
        done

        trimmed_alignment_size_gappyout=$(wc --line < SoP_by_column_${filename}_trimmed_gappyout.txt)
        if [[ $(echo "$trimmed_alignment_size_gappyout % 2" | bc) == '0' ]]; then
            median_pos_first=$(echo "$trimmed_alignment_size_gappyout / 2" | bc)
            median_pos_second=$(echo "$median_pos_first + 1" | bc)
            median_first=$(sed "${median_pos_first}q;d" SoP_by_column_${filename}_trimmed_gappyout.txt)
            median_second=$(sed "${median_pos_second}q;d" SoP_by_column_${filename}_trimmed_gappyout.txt)
            median_alignment_gappyout=$(echo "($median_first + $median_second) / 2" | bc)
        else
            median_pos=$(echo "$trimmed_alignment_size_gappyout / 2 + 1" | bc)
            median_alignment_gappyout=$(sed "${median_pos}q;d" SoP_by_column_${filename}_trimmed_gappyout.txt)
        fi

        min_SoP_trimmed_alignment_gappyout=$(sort SoP_by_column_${filename}_trimmed_gappyout.txt | head -n 1)
        max_SoP_trimmed_alignment_gappyout=$(sort SoP_by_column_${filename}_trimmed_gappyout.txt | tail -n 1)
        SoP_by_column_trimmed_sum=$(paste -sd+ SoP_by_column_${filename}_trimmed_gappyout.txt | bc)
        trimmed_alignment_SoP_gappyout=$(echo "$SoP_by_column_trimmed_sum / $trimmed_alignment_size_gappyout" | bc -l)
        std_alignment_SoP_gappyout=$(awk -v avg_SoP=$trimmed_alignment_SoP_gappyout -v size=$trimmed_alignment_size_gappyout '{n += ($1 - avg_SoP) ** 2}; END{print sqrt(n/size)}'\
        SoP_by_column_${filename}_trimmed_gappyout.txt)
        rm SoP_by_column_${filename}_trimmed_gappyout.txt 
    else
        trimmed_alignment_size_gappyout=0
        median_alignment_gappyout=0
        min_SoP_trimmed_alignment_gappyout=0
        max_SoP_trimmed_alignment_gappyout=0
        std_alignment_SoP_gappyout=0
        trimmed_alignment_SoP_gappyout=0
    fi



    trimmed_column_positions_strictplus=$($trimal_local -in ${file}_no_all_gaps_cols -strictplus -colnumbering | grep "#ColumnsMap" | awk '{print substr($0,13)}')
    if [ -n "$trimmed_column_positions_strictplus" ]; then
        > SoP_by_column_${filename}_trimmed_strictplus.txt
        # -colnumbering returns columns starting by 0 and -sfc starting by 1
        for column_position in ${trimmed_column_positions_strictplus//,/};
        do
            col_pos=$((column_position + 1))
            sed "${col_pos}q;d" SoP_by_column_$filename.txt >> SoP_by_column_${filename}_trimmed_strictplus.txt
        done

        trimmed_alignment_size_strictplus=$(wc --line < SoP_by_column_${filename}_trimmed_strictplus.txt)
        if [[ $(echo "$trimmed_alignment_size_strictplus % 2" | bc) == '0' ]]; then
            median_pos_first=$(echo "$trimmed_alignment_size_strictplus / 2" | bc)
            median_pos_second=$(echo "$median_pos_first + 1" | bc)
            median_first=$(sed "${median_pos_first}q;d" SoP_by_column_${filename}_trimmed_strictplus.txt)
            median_second=$(sed "${median_pos_second}q;d" SoP_by_column_${filename}_trimmed_strictplus.txt)
            median_alignment_strictplus=$(echo "($median_first + $median_second) / 2" | bc)
        else
            median_pos=$(echo "$trimmed_alignment_size_strictplus / 2 + 1" | bc)
            median_alignment_strictplus=$(sed "${median_pos}q;d" SoP_by_column_${filename}_trimmed_strictplus.txt)
        fi

        min_SoP_trimmed_alignment_strictplus=$(sort SoP_by_column_${filename}_trimmed_strictplus.txt | head -n 1)
        max_SoP_trimmed_alignment_strictplus=$(sort SoP_by_column_${filename}_trimmed_strictplus.txt | tail -n 1)
        SoP_by_column_trimmed_sum=$(paste -sd+ SoP_by_column_${filename}_trimmed_strictplus.txt | bc)
        trimmed_alignment_SoP_strictplus=$(echo "$SoP_by_column_trimmed_sum / $trimmed_alignment_size_strictplus" | bc -l)
        std_alignment_SoP_strictplus=$(awk -v avg_SoP=$trimmed_alignment_SoP_strictplus -v size=$trimmed_alignment_size_strictplus '{n += ($1 - avg_SoP) ** 2}; END{print sqrt(n/size)}'\
        SoP_by_column_${filename}_trimmed_strictplus.txt)
        rm SoP_by_column_${filename}_trimmed_strictplus.txt 
    else
        trimmed_alignment_size_strictplus=0
        median_alignment_strictplus=0
        min_SoP_trimmed_alignment_strictplus=0
        max_SoP_trimmed_alignment_strictplus=0
        std_alignment_SoP_strictplus=0
        trimmed_alignment_SoP_strictplus=0
    fi

    rm ${file}_filtered_refs
    rm ${file}_no_all_gaps_cols
    rm SoP_output_${filename}.txt
    rm SoP_by_column_${filename}.txt

    echo "$filename,$original_num_seqs,$original_num_cols,$clean_num_seqs,$clean_num_cols,$min_SoP_alignment,"\
    "$max_SoP_alignment,$median_alignment,$std_alignment,$alignment_SoP,$trimmed_alignment_size_gappyout,"\
    "$min_SoP_trimmed_alignment_gappyout,$max_SoP_trimmed_alignment_gappyout,$median_alignment_gappyout,"\
    "$std_alignment_SoP_gappyout,$trimmed_alignment_SoP_gappyout,$trimmed_alignment_size_strictplus,"\
    "$min_SoP_trimmed_alignment_strictplus,$max_SoP_trimmed_alignment_strictplus,$median_alignment_strictplus,"\
    "$std_alignment_SoP_strictplus,$trimmed_alignment_SoP_strictplus"
}

compute_sum_of_pairs() {
    file=$1
    filename=$(echo $file | awk -F '/' '{print $NF}')
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    #echo "trimming $filename"

    # write ref alignment in set_alignments_SoP.txt
    echo $ref_alignment_name.ref > ref_alignment_SoP_$ref.txt
    $trimal_local -sfc -compareset ref_alignment_SoP_$ref.txt -forceselect $file > SoP_output_$filename.txt
    original_num_cols=$(cat SoP_output_$filename.txt | tail -n 1 | awk '{print $1}')
    tail -n +9 SoP_output_${filename}.txt | awk '{print $2}' > SoP_by_column_${filename}.txt
    SoP_by_column_sum=$(paste -sd+ SoP_by_column_${filename}.txt | bc)
    clean_alignment_size=$(wc --line < SoP_by_column_${filename}.txt)
    alignment_SoP=$(echo "$SoP_by_column_sum / $clean_alignment_size" | bc -l)

    trimmed_column_positions_gappyout=$($trimal_local -in $file -gappyout -colnumbering | grep "#ColumnsMap" | awk '{print substr($0,13)}')
    > SoP_by_column_${filename}_trimmed_gappyout.txt
    # -colnumbering returns columns starting by 0 and -sfc starting by 1
    for column_position in ${trimmed_column_positions_gappyout//,/};
    do
        col_pos=$((column_position + 1))
        sed "${col_pos}q;d" SoP_by_column_$filename.txt >> SoP_by_column_${filename}_trimmed_gappyout.txt
    done


    SoP_by_column_trimmed_sum=$(paste -sd+ SoP_by_column_${filename}_trimmed_gappyout.txt | bc)
    trimmed_alignment_size_gappyout=$(wc --line < SoP_by_column_${filename}_trimmed_gappyout.txt)
    trimmed_alignment_SoP_gappyout=$(echo "$SoP_by_column_trimmed_sum / $trimmed_alignment_size_gappyout" | bc -l)

    trimmed_column_positions_strictplus=$($trimal_local -in $file -strictplus -colnumbering | grep "#ColumnsMap" | awk '{print substr($0,13)}')
    > SoP_by_column_${filename}_trimmed_strictplus.txt
    # -colnumbering returns columns starting by 0 and -sfc starting by 1
    for column_position in ${trimmed_column_positions_strictplus//,/};
    do
        col_pos=$((column_position + 1))
        sed "${col_pos}q;d" SoP_by_column_$filename.txt >> SoP_by_column_${filename}_trimmed_strictplus.txt
    done


    SoP_by_column_trimmed_sum=$(paste -sd+ SoP_by_column_${filename}_trimmed_strictplus.txt | bc)
    trimmed_alignment_size_strictplus=$(wc --line < SoP_by_column_${filename}_trimmed_strictplus.txt)
    trimmed_alignment_SoP_strictplus=$(echo "$SoP_by_column_trimmed_sum / $trimmed_alignment_size_strictplus" | bc -l)

    echo "$filename; clean_alignment_columns=$clean_alignment_size; clean_msa_alignment_SoP=$alignment_SoP; trimmed_alignment_SoP_gappyout=$trimmed_alignment_SoP_gappyout;"\
    "trimmed_alignment_SoP_strictplus=$trimmed_alignment_SoP_strictplus"

    echo $filename >> alignments.txt
}

remove_all_gaps_columns() {
    file=$1
    $trimal_local -in $file -noallgaps > "${file}_no_all_gaps_cols"
    echo "Processed $file"
}

filter_reference_sequences() {
    file=$1
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    python3 $filter_reference_sequences_script -i $file -r $ref_alignment_name.ref
    echo "Processed $file"
}

compute_sum_of_pairs_stats $file
