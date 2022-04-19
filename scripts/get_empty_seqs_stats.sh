#!/bin/bash

start_time=$(date +%s.%N)
empty_seqs_statistics="empty_seqs_statistics"
if [ ! -d "$empty_seqs_statistics" ]; then
    mkdir $empty_seqs_statistics
fi
cd $empty_seqs_statistics

> residue_type.txt
> taxon.txt
> problem_num.txt
> error_problem.txt
> MSA_tools.txt
> MSA_filter_tools.txt
> log.txt
> log_writing.txt


write_results() {
    file="../all_gaps_indets_alignments.txt"
    while read -r filename; do
        filename=$(echo $filename | awk -F '/' '{print $(NF-3)"_"$(NF-2)"_"$NF}')
		echo "writing $filename" >> log.txt &
		if [[ $filename == *.fa_err ]]; then
			echo "True" >> error_problem.txt &
		else
			echo "False" >> error_problem.txt &
		fi
		echo $filename | awk -F '_' '{print $1}' >> residue_type.txt &
		echo $filename | awk -F '_' '{print $2}' >> taxon.txt &
		echo $filename | awk -F '.' '{print $1}' | grep -o '[0-9]\+' >> problem_num.txt &
		if [[ $filename != *".seqs."* ]]; then
			MSA_tool=$(echo $filename | awk -F '.' '{print $2}')
			MSA_filter_tool=""
			if [[ $MSA_tool == Guidance* ]]; then
				MSA_tool=$(echo $MSA_tool | awk -F 'Guidance' '{print $2}')
				echo $MSA_tool  >> MSA_tools.txt &
				MSA_filter_tool='Guidance'
				echo $MSA_filter_tool >> MSA_filter_tools.txt &
			else
				echo $MSA_tool >> MSA_tools.txt &
				MSA_filter_tool=$(echo $filename | awk -F '.' '{if (NF == 4) print $3; else print "None";}')
				echo $MSA_filter_tool >> MSA_filter_tools.txt &
			fi

		else
			echo "Original" >> MSA_tools.txt &
			echo "Original" >> MSA_filter_tools.txt &
		fi
		echo "Finished with $filename" >> log_writing.txt &
		wait
    done < $file
}


write_results

end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff >> log.txt
