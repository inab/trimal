#!/bin/bash

seed_alignments=$(find trimal/dataset/*.fasta -not -name '*014*' -not -name '*015*' -not -name '*016*' -not -name '*018*' -not -name '*019*' -not -name '*020*' -not -name '*021*' -not -name '*023*' -not -name '*025*' -not -name '*026*' -not -name '*027*' -not -name '*028*' -not -name '*029*' -not -name '*030*' -not -name '*031*' -not -name '*032*' -not -name '*033*' -not -name '*034*' -not -name '*035*' -not -name '*037*' -not -name '*039*' -not -name '*040*' -not -name '*051*' -not -name '*052*' -not -name '*054*' -not -name '*057*' -not -name '*060*'  -not -name '*062*' -not -name '*064*' -not -name '*067*' -not -name '*070*' -not -name '*073*' -not -name '*076*' -not -name '*088*' -not -name '*090*' -not -name '*091*' -not -name '*092*' -not -name '*093*')

seed_alignments_count=$(find trimal/dataset/*.fasta -not -name '*014*' -not -name '*015*' -not -name '*016*' -not -name '*018*' -not -name '*019*' -not -name '*020*' -not -name '*021*' -not -name '*023*' -not -name '*025*' -not -name '*026*' -not -name '*027*' -not -name '*028*' -not -name '*029*' -not -name '*030*' -not -name '*031*' -not -name '*032*' -not -name '*033*' -not -name '*034*' -not -name '*035*' -not -name '*037*' -not -name '*039*' -not -name '*040*' -not -name '*051*' -not -name '*052*' -not -name '*054*' -not -name '*057*' -not -name '*060*'  -not -name '*062*' -not -name '*064*' -not -name '*067*' -not -name '*070*' -not -name '*073*' -not -name '*076*' -not -name '*088*' -not -name '*090*' -not -name '*091*' -not -name '*092*' -not -name '*093*' | wc -l)


for alignment in $seed_alignments
do
	seqs=100
	while [ $seqs -lt 100000 ]
	do
		cols=100
		while [ $cols -lt 1000000 ]
		do
			alignment_name=$(basename $alignment .fasta)
			echo "python3 ../trimal/scripts/generateRandomAlignmentsUsingAsSeedRealAlignments.py -i ../$alignment -s $seqs -r $cols -m 1000 > simulated_alignments/${alignment_name}.${seqs}.${cols}.fasta"
			((cols*=2))
		done
		((seqs*=2))
	done
done
