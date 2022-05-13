#!/bin/bash

file=$1

dataset_ref_alignments="cedric_datasets/refs"
trimal_local="trimal/bin/trimal"
filter_reference_sequences_script="trimal/scripts/filter_homfam_alignments.py"


compute_sum_of_pairs_stats() {
    file=$1
    original_num_seqs=$(grep -c ">" $file)
    original_num_cols=$(while read line; do echo $line; [ -z $line ] && break; done <<<  "$(awk -F '>' '{print $1}' ${file} | tail -n +2)" | paste -sd '' | wc -c)
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    python3 $filter_reference_sequences_script -i $file -r $ref_alignment_name.ref
    $trimal_local -in ${file}_filtered_refs -noallgaps > "${file}_no_all_gaps_cols"
    clean_num_seqs=$(grep -c ">" ${file}_no_all_gaps_cols)
    clean_num_cols=$(while read line; do echo $line; [ -z $line ] && break; done <<<  "$(awk -F '>' '{print $1}' ${file}_no_all_gaps_cols | tail -n +2)" | paste -sd '' | wc -c)

    filename=$(echo $file | awk -F '/' '{print $NF}')
    ref=$(basename $file | awk -F '.' '{print $1}')
    ref_alignment_name="$dataset_ref_alignments/$ref"
    if [ ! -s ref_alignment_SoP_$ref.txt ]; then
        echo $ref_alignment_name.ref > ref_alignment_SoP_$ref.txt
    fi
    
    $trimal_local -sfc -compareset ref_alignment_SoP_$ref.txt -forceselect ${file}_no_all_gaps_cols > SoP_output_${filename}.txt
    tail -n +9 SoP_output_${filename}.txt | awk '{print $2}' > SoP_by_column_${filename}.txt
    SoP_by_column_sum=$(paste -sd+ SoP_by_column_${filename}.txt | bc)
    alignment_SoP=$(echo "$SoP_by_column_sum / $clean_num_cols" | bc -l)


    trimmed_column_positions_gappyout=$($trimal_local -in ${file}_no_all_gaps_cols -gappyout -colnumbering | grep "#ColumnsMap" | awk '{print substr($0,13)}')
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


    trimmed_column_positions_strictplus=$($trimal_local -in ${file}_no_all_gaps_cols -strictplus -colnumbering | grep "#ColumnsMap" | awk '{print substr($0,13)}')
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

    rm ${file}_filtered_refs
    rm ${file}_no_all_gaps_cols
    rm SoP_output_${filename}.txt
    rm SoP_by_column_${filename}.txt
    rm SoP_by_column_${filename}_trimmed_gappyout.txt 
    rm SoP_by_column_${filename}_trimmed_strictplus.txt

    echo "$filename,$original_num_seqs,$original_num_cols,$clean_num_seqs,$clean_num_cols,$alignment_SoP,"\
    "$trimmed_alignment_size_gappyout,$trimmed_alignment_SoP_gappyout,$trimmed_alignment_size_strictplus,"\
    "$trimmed_alignment_SoP_strictplus"
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
