#!/bin/bash

start_time=$(date +%s.%N)
dataset_alignments="../cedric_datasets/homfam/alignments"
dataset_ref_alignments="../cedric_datasets/homfam/refs"
homfam_filtered="homfam_filtered"
trimal_local="../trimal/bin/trimal"
if [ ! -d "$homfam_filtered" ]; then
    mkdir $homfam_filtered
fi
cd $homfam_filtered
> log.txt
> sum_of_pairs.txt
> alignments.txt

trim_alignments() {
    file=$1
    filename=$(echo $file | awk -F '/' '{print $NF}')
    ref_alignment_name=$(echo $filename | awk -F '.' '{print $1}')
    echo "trimming $filename" >> log.txt
    $trimal_local -in $file -gappyout > gappyout_$filename &
    $trimal_local -in $file -strict > strict_$filename &
    $trimal_local -in $file -strictplus > strictplus_$filename &
    $trimal_local -in $file -automated1 > automated1_$filename &
    wait
    (t_coffee -other_pg aln_compare -al1 $dataset_ref_alignments/$ref_alignment_name.ref -al2 $file -compare_mode sp \
        | grep -v "seq1" | grep -v '*' | awk '{print $4}' ORS="\t" >> sum_of_pairs.txt)
    echo $filename >> alignments.txt
    (t_coffee -other_pg aln_compare -al1 $dataset_ref_alignments/$ref_alignment_name.ref -al2 gappyout_$filename -compare_mode sp \
        | grep -v "seq1" | grep -v '*' | awk '{print $4}' ORS="\t" >> sum_of_pairs.txt)
    echo "gappyout_$filename" >> alignments.txt
    (t_coffee -other_pg aln_compare -al1 $dataset_ref_alignments/$ref_alignment_name.ref -al2 strict_$filename -compare_mode sp \
        | grep -v "seq1" | grep -v '*' | awk '{print $4}' ORS="\t" >> sum_of_pairs.txt)
    echo "strict_$filename" >> alignments.txt
    (t_coffee -other_pg aln_compare -al1 $dataset_ref_alignments/$ref_alignment_name.ref -al2 strictplus_$filename -compare_mode sp \
        | grep -v "seq1" | grep -v '*' | awk '{print $4}' ORS="\t" >> sum_of_pairs.txt)
    echo "strictplus_$filename" >> alignments.txt
    (t_coffee -other_pg aln_compare -al1 $dataset_ref_alignments/$ref_alignment_name.ref -al2 automated1_$filename -compare_mode sp \
        | grep -v "seq1" | grep -v '*' | awk '{print $4}' ORS="\t" >> sum_of_pairs.txt)
    echo "automated1_$filename" >> alignments.txt
    echo "trimmed $filename" >> log.txt
}

counter=0
for file in $dataset_alignments/*
do
    trim_alignments $file &
    ((counter++))
    if [[ $counter -eq 5 ]]; then
        wait
        counter=0
        break
    fi
done



end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
