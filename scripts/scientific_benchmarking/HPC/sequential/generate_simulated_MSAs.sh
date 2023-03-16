#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=48
#SBATCH --qos=debug
#SBATCH --output=generate_simulated_MSAs.txt
#SBATCH --error=generate_simulated_MSAs.err
module load python

seqs=100
counter=1
while [ $seqs -lt 1000 ]
do
	cols=100
	while [ $cols -lt 1000 ]
	do
		alignment=$(ls "cedric_datasets/alignments" | sed $counter"q;d")
		alignment_name=$(basename $alignment .fasta)
		python3 trimal/scripts/generateRandomAlignmentsUsingAsSeedRealAlignments.py -i cedric_datasets/alignments/$alignment -s $seqs -r $cols -m 1000 > simulated_alignments/${alignment_name}.${seqs}.${cols}.fasta &
		echo "Generating ${alignment_name}.${seqs}.${cols}.fasta"
		((counter+=5))
		((cols*=2))
	done
	((seqs*=2))
done

wait
