#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=24
#SBATCH --array=401-500
#SBATCH --qos=bsc_ls
#SBATCH --output=automated2_arrays.txt
module load iq-tree

start_time=$(date +%s.%N)
dataset="dessimoz/02.Data/SpeciesTreeDiscordanceTest.Enriched"
problems_to_ignore="all_gaps_indets_alignments_unique_AA.txt"

counter=0
id=$(printf %03d $SLURM_ARRAY_TASK_ID)
for problem in $dataset/AA/Bacteria/problem${id}*
do
	problem_name=$(basename $problem | awk -F '_' '{print $1}')
	residue_taxon_problem=$(echo $problem | awk -F '/' '{print $(NF-2)"_"$(NF-1)"_"$NF}')
	if grep -Fxq ${residue_taxon_problem} ${problems_to_ignore}; then
		echo "Ignored $residue_taxon_problem"
	else
		for file in $problem/*{ClustalW2,Mafft,Prank,T-Coffee}*trimAl_automated2*.fa
		do
			echo "Calculating tree for ${file}"
			iqtree -nt 2 -quiet -redo -mredo -cptime 240 -mem 4G -cmin 4 -cmax 10 -s $file -bb 1000 -mset WAG,LG,JTT &
			((counter++))
			if [ $counter -eq 24 ]; then
				wait
				counter=0
			fi
		done
	fi
done

end_time=$(date +%s.%N)
diff=$(echo "$end_time - $start_time" | bc)
echo $diff
