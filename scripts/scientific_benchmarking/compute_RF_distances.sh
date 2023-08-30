#!/bin/bash

file=$1

task() {
    file=$1
    file_basename=$(basename $file)
    start_file_time=$(date +%s.%N)
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    folder_path=$(dirname $file)
    tree_file="$file.treefile"
    #if [ ! -s "$tree_file" ]; then
    #    iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $file -bb 1000 -mset GTR
    #fi

    number_sequences=$(grep -c ">" $file)
    number_columns=$(awk '/^>/{l=0; next}{l+=length($0)}END{print l}' $file)
    ref_tree=$(find $folder_path -name *.reference.nwk)
    RF_distance=$(ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename.txt &&
            awk 'FNR == 3 {print $9}' temp_RF_distance_$filename.txt || 
            echo -1)
    rm temp_RF_distance_$filename.txt
    residue_type=$(echo $filename | awk -F '_' '{print $1}')
	taxon=$(echo $filename | awk -F '_' '{print $2}')
	problem_num=$(echo $filename | awk -F '.' '{print $1}' | grep -o '[0-9]\+')

    if [[ $filename != *".seqs."* && $filename != *.fa_err  ]]; then
        MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
        MSA_filter_tool=""
        if [[ $MSA_tool == Guidance* ]]; then
            MSA_tool=$(echo $MSA_tool | awk -F 'Guidance' '{print $2}')
            MSA_filter_tool='Guidance'
        else
            MSA_filter_tool=$(echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}')
            MSA_tool_filename=$(echo $file_basename | awk -F '.' '{print $1"."$2".fa"}')
        fi
        RF_distance_diff=-1
        if [[ $MSA_filter_tool != "None" ]]; then
            if [[ $MSA_filter_tool == "Guidance" ]]; then
                if [[ $MSA_tool == "ClustalW" ]]; then
                    problem_name=$(echo $file_basename | awk -F '.' '{print $1}')
                    MSA_tool_filename="$problem_name.${MSA_tool}2.fa"
                else
                    problem_name=$(echo $file_basename | awk -F '.' '{print $1}')
                    MSA_tool_filename="$problem_name.${MSA_tool}.fa"
                fi
            else
                file_suffix=$(echo $file_basename | awk -F '.' '{print $(NF-3)"."$(NF-2)"."$NF}')
                MSA_tool_filename="$file_suffix"
            fi

            tree_file="$folder_path/$MSA_tool_filename.treefile"

            if [ -s RF_distances_temp_${taxon}/temp_RF_distance_$MSA_tool_filename.txt ]; then
                RF_distance_tool=$(awk 'FNR == 3 {print $9}' RF_distances_temp_${taxon}/temp_RF_distance_$MSA_tool_filename.txt)
            else
                RF_distance_tool=$(ete3 compare -t $tree_file -r $ref_tree --unrooted > RF_distances_temp_${taxon}/temp_RF_distance_$MSA_tool_filename.txt &&
                awk 'FNR == 3 {print $9}' RF_distances_temp_${taxon}/temp_RF_distance_$MSA_tool_filename.txt || 
                echo -1)
            fi

            RF_distance_diff=$(echo "$RF_distance_tool - $RF_distance" | bc)
        fi
        if [[ -z $RF_distance_diff ]]; then
            RF_distance_diff=-100
        fi
    else
        MSA_tool="Original"
        MSA_filter_tool="Original"
        number_columns=-1
        RF_distance=-1
        RF_distance_diff=-1
    fi

    if (( $(echo "$RF_distance_diff == 0" | bc -l) )) ; then
        RF_distance_change="unchanged"
    elif (( $(echo "$RF_distance_diff < 0" | bc -l) )); then
        RF_distance_change="worse"
    elif (( $(echo "$RF_distance_diff > 0" | bc -l) )); then
        RF_distance_change="better"
    else
        RF_distance_change="error"
    fi

    echo "$taxon,$problem_num,$number_sequences,$number_columns,$MSA_tool,$MSA_filter_tool,$RF_distance,$RF_distance_diff,$RF_distance_change"

    #rm temp_RF_distance_$MSA_tool_filename.txt
}


task "$file"