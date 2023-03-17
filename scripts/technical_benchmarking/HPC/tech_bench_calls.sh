#!/bin/bash

dataset="simulated_alignments/chunk1/from_5GB_to_20GB"

for alignment in $dataset/*
do
	for rep in {1..100}
	do
		echo "./../trimal/scripts/technical_performance_tests.sh $alignment $rep"
	done
done
