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
> avg_gaps.txt
> avg_seq_identity.txt
> RF_distance.txt
> residue_type.txt
> taxon.txt
> problem_num.txt
> MSA_tools.txt
> MSA_filter_tools.txt
> log.txt

task(){
    start_file_time=$(date +%s.%N)
    file=$1
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    folder_path=$(dirname $file)
    filename_no_ext=$(basename $file .fa)
    tree_file="$folder_path/$filename_no_ext.nwk"
    if [[ $filename != *".seqs."* ]]; then
        $trimal_local -in $file -sgc > temporal_out_$filename.txt
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
        $trimal_local -in $file -sident | grep "## AverageIdentity" | awk '{print $3}' > avg_seq_identity_$filename.txt &
        if [ -s $tree_file ]; then
            ref_tree=$(find $folder_path -name *.reference.nwk)
            ete3 compare -t $tree_file -r $ref_tree | awk 'FNR == 3 {print $9}' > RF_distance_$filename.txt &
        fi
    fi
    wait
    end_file_time=$(date +%s.%N)
    file_time=$(echo "$end_file_time - $start_file_time" | bc)
    filename="${file##*/}"
    #echo "$filename took $file_time  seconds"
}


N=5
for residue_type in $dataset/*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/*
        do
            for problem in $taxon/problem0001*
            do
                num_seq=$(grep ">" "$problem/problem0001.seqs.fa" | wc -l)
                for file in $problem/*.fa
                do
                    ((i=i%N)); ((i++==0)) && wait
                    echo $num_seq >> number_sequences.txt
                    task "$file"  &
                done
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
            for problem in $taxon/problem0001*
            do
                for file in $problem/*.fa
                do
                    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
                    if [[ $filename != *".seqs."* ]]; then
                        echo $residue_type | awk -F '/' '{print $(NF)}' >> residue_type.txt
                        echo $taxon | awk -F '/' '{print $(NF)}' >> taxon.txt 
                        echo $problem | awk -F '/' '{print $(NF)}' >> problem_num.txt
                        MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
                        if  [[ $MSA_tool == Guidance* ]]; then
                            echo $MSA_tool | awk -F 'Guidance' '{print $2}' >> MSA_tools.txt
                            echo 'Guidance' >> MSA_filter_tools.txt
                        else
                            echo $MSA_tool >> MSA_tools.txt
                            echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}' >> MSA_filter_tools.txt
                        fi
                        cat blocks_outputs_$filename.txt >> blocks_outputs.txt && rm blocks_outputs_$filename.txt
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
                        echo $residue_type | awk -F '/' '{print $(NF)}' >> residue_type.txt
                        echo $taxon | awk -F '/' '{print $(NF)}' >> taxon.txt 
                        echo $problem | awk -F '/' '{print $(NF)}' >> problem_num.txt
                        echo "Original" >> MSA_tools.txt
                        echo "Original" >> MSA_filter_tools.txt
                        echo -1 >> blocks_outputs.txt
                        echo -1 >> number_columns.txt
                        echo -1 >> avg_gaps.txt
                        echo -1 >> avg_seq_identity.txt
                        echo -1 >> RF_distance.txt
                    fi
                    echo "Finished with $filename" >> log.txt
                done
            done
        done
    fi
done


end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
echo $diff