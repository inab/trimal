#!/bin/bash

res_type=$1
problem_nums=$2
taxon_name=$3
problems_to_ignore="../all_gaps_indets_alignments_unique_${res_type}_${taxon_name}_${problem_nums}.txt"

start_time=$(date +%s.%N)
scripts="../trimal/scripts"
dataset="../dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
alignment_statistics="automated2_trees"
trimal_local="../trimal/bin/trimal"
#problems_to_ignore="../all_gaps_indets_alignments_unique_$res_type.txt"
if [ ! -d "$alignment_statistics" ]; then
    mkdir $alignment_statistics
fi
cd $alignment_statistics

> residue_type.txt
> taxon.txt
> problem_num.txt
> MSA_tools.txt
> MSA_filter_tools.txt
> RF_distance.txt
> RF_distance_diff.txt
> log.txt
> filenames_to_write.txt


trim_alignment(){
	input_file=$1
	method=$2
	min_percentage_columns=$3
	min_columns=$4
	output_file_prefix=$5
	min_percentage_columns_num=$(echo 0.$min_percentage_columns | bc)

	(python3 $scripts/automated2_algorithm.py -i $input_file -a $method --min_percentage_columns $min_percentage_columns_num \
		--min_columns $min_columns > $output_file_prefix.trimAl_automated2_${method}_${min_percentage_columns}_${min_columns}.fa)
}

trim_alignments(){
	file=$1
    folder_path=$(dirname $file)
	file_basename=$(basename $file .fa)

	echo "Calculating trimmed alignment for ${file_basename}" >> log.txt

	min_percentage_columns_list=(30 35 40 45)
	min_columns_list=(100 150 200 250)
	for min_percentage_columns in ${min_percentage_columns_list[*]}; do
		for min_columns in ${min_columns_list[*]}; do
			trim_alignment $file "gaps" $min_percentage_columns $min_columns "$folder_path/${file_basename}" &
			trim_alignment $file "similarity" $min_percentage_columns $min_columns "$folder_path/${file_basename}" &
			trim_alignment $file "combined" $min_percentage_columns $min_columns "$folder_path/${file_basename}" &
		done
	done
}

generate_tree(){
	output_file_prefix=$1
	method=$2
	min_percentage_columns=$3
	min_columns=$4
	iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $output_file_prefix.trimAl_automated2_$method_${min_percentage_columns}_${min_columns}.fa -bb 1000 -mset WAG,LG,JTT
}

generate_trees(){
	file=$1
    folder_path=$(dirname $file)
	file_basename=$(basename $file .fa)

	for min_percentage_columns in ${min_percentage_columns_list[*]}; do
		for min_columns in ${min_columns_list[*]}; do
			generate_tree "$folder_path/${file_basename}" "gaps" $min_percentage_columns $min_columns &
			generate_tree "$folder_path/${file_basename}" "similarity" $min_percentage_columns $min_columns &
			generate_tree "$folder_path/${file_basename}" "combined" $min_percentage_columns $min_columns &
		done
	done
}

task(){
    file=$1
    start_file_time=$(date +%s.%N)
    filename=$(echo $file | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
    folder_path=$(dirname $file)
	file_basename=$(basename $file .fa)

	echo "Calculating trimmed alignment for ${file_basename}" >> log.txt

	# generate trimmed alignments


	#(python3 $scripts/automated2_algorithm.py -i $file -a gaps --min_percentage_columns 0.3 --min_columns 100 > $folder_path/${file_basename}.trimAl_automated2_gaps_30_100.fa) &

	#(python3 $scripts/automated2_algorithm.py -i $file -a similarity > $folder_path/${file_basename}.trimAl_automated2_sim.fa) &

	#(python3 $scripts/automated2_algorithm.py -i $file -a combined > $folder_path/${file_basename}.trimAl_automated2_comb.fa) &
	wait
	return

	echo "Calculating tree for ${file_basename}.trimAl_automated2.fa" >> log.txt
	# generate new trees


	#iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $folder_path/${file_basename}.trimAl_automated2_gaps.fa -bb 1000 -mset WAG,LG,JTT &
	#iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $folder_path/${file_basename}.trimAl_automated2_sim.fa -bb 1000 -mset WAG,LG,JTT &
	#iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $folder_path/${file_basename}.trimAl_automated2_comb.fa -bb 1000 -mset WAG,LG,JTT &
	wait

	return
	
	tree_file="$file.treefile"
	ref_tree=$(find $folder_path -name *.reference.nwk)
    (ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename.txt &&
    awk 'FNR == 3 {print $9}' temp_RF_distance_$filename.txt || 
    echo -1) > RF_distance_$filename.txt
	rm temp_RF_distance_$filename.txt

	filename_automated2_gaps=$(basename $filename .fa)
	filename_automated2_gaps=${filename_automated2_gaps}.trimAl_automated2_gaps.fa
	tree_file="$folder_path/${file_basename}.trimAl_automated2_gaps.fa.treefile"
    (ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename_automated2_gaps.txt &&
    awk 'FNR == 3 {print $9}' temp_RF_distance_$filename_automated2_gaps.txt || 
    echo -1) > RF_distance_$filename_automated2_gaps.txt
	rm temp_RF_distance_$filename_automated2_gaps.txt
	echo $filename_automated2_gaps >> filenames_to_write.txt

	filename_automated2_sim=$(basename $filename .fa)
	filename_automated2_sim=${filename_automated2_sim}.trimAl_automated2_sim.fa
	tree_file="$folder_path/${file_basename}.trimAl_automated2_sim.fa.treefile"
    (ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename_automated2_sim.txt &&
    awk 'FNR == 3 {print $9}' temp_RF_distance_$filename_automated2_sim.txt || 
    echo -1) > RF_distance_$filename_automated2_sim.txt
	rm temp_RF_distance_$filename_automated2_sim.txt
	echo $filename_automated2_sim >> filenames_to_write.txt

	filename_automated2_comb=$(basename $filename .fa)
	filename_automated2_comb=${filename_automated2_comb}.trimAl_automated2_comb.fa
	tree_file="$folder_path/${file_basename}.trimAl_automated2_comb.fa.treefile"
    (ete3 compare -t $tree_file -r $ref_tree --unrooted > temp_RF_distance_$filename_automated2_comb.txt &&
    awk 'FNR == 3 {print $9}' temp_RF_distance_$filename_automated2_comb.txt || 
    echo -1) > RF_distance_$filename_automated2_comb.txt
	rm temp_RF_distance_$filename_automated2_comb.txt
	echo $filename_automated2_comb >> filenames_to_write.txt
    
    
    end_file_time=$(date +%s.%N)
    file_time=$(echo "$end_file_time - $start_file_time" | bc)
    echo "$filename processed in $file_time seconds"
    echo "$filename processed in $file_time seconds" >> log.txt 
}


write_results() {
    file="filenames_to_write.txt"
    while read -r filename; do
		echo "writing $filename" >> log.txt &
		echo $filename | awk -F '_' '{print $1}' >> residue_type.txt &
		echo $filename | awk -F '_' '{print $2}' >> taxon.txt &
		echo $filename | awk -F '.' '{print $1}' | grep -o '[0-9]\+' >> problem_num.txt &
		MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
		echo $MSA_tool >> MSA_tools.txt &
		echo "trimAl_automated2" >> MSA_filter_tools.txt &

		RF_distance_diff=-1

		file_suffix=$(echo $filename | awk -F '.' '{print $(NF-3)"."$(NF-2)"."$NF}')

		RF_distance_MSA_tool_file="RF_distance_$file_suffix.txt"
		RF_distance_filter_tool=$(cat RF_distance_$filename.txt)
		RF_distance_tool=$(cat $RF_distance_MSA_tool_file)
		RF_distance_diff=$(echo "$RF_distance_tool - $RF_distance_filter_tool" | bc)

		(cat RF_distance_$filename.txt || echo -1) >> RF_distance.txt

		if [[ -z $RF_distance_diff ]]; then
			(echo -100 >> RF_distance_diff.txt) &
		else
			(echo $RF_distance_diff >> RF_distance_diff.txt) &
		fi
		
		wait
    done < $file

	if [ -s filenames_to_write.txt ]; then
		rm RF_distance_*problem*.txt
    fi

    > filenames_to_write.txt
}


counter=0
line_number=1
for residue_type in $dataset/$res_type*
do
    if [ -d $residue_type ]; then
        for taxon in $residue_type/$taxon_name*
        do
	    	residue_taxon_problem=""
            for problem in $taxon/problem$problem_nums*
            do
                problem_name=$(basename $problem | awk -F '_' '{print $1}')
				residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
				residue_taxon_problem_to_ignore=$(sed "${line_number}q;d" $problems_to_ignore)
				if [[ $residue_taxon_problem != $residue_taxon_problem_to_ignore ]]; then
					for file in $problem/*.fa
					do
						if [[ $file == *"${problem_name}.ClustalW2.fa"* || $file == *"${problem_name}.Mafft.fa"* || $file == *"${problem_name}.Prank.fa"* || $file == *"${problem_name}.T-Coffee.fa"* ]]; then
							#trim_alignments "$file"  &
							generate_trees "$file"  &
						fi
					done
				else
					((line_number++))
					echo "Ignored $residue_taxon_problem" >> log.txt
				fi

				((counter++))
				if [ $counter -eq 2 ]; then
					wait
					counter=0
					#write_results
				fi
            done
		done
    fi
done

wait
#write_results

end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
