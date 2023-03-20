#!/bin/bash

# Generate the commands to execute
dataset="02.Data/SpeciesTreeDiscordanceTest.Enriched"
number_commands=0
> commands_to_execute.txt

for residue_type in $dataset/{AA,DNA}
do
	for species_type in $residue_type/*
	do
		for problem in $species_type/problem*
		do
			echo "command $problem" >> commands_to_execute.txt
			((number_commands++))
		done
	done
done

# Split commands into several files
number_jobs=300
commands_per_job=$((number_commands / number_jobs))
mkdir commands_to_execute
split --lines $commands_per_job commands_to_execute.txt commands_to_execute_split/commands_to_execute_

# Generate the greasy commands
> commands_to_execute_greasy.txt
for file in commands_to_execute_split/*
do
	echo "/apps/GREASY/2.2/INTEL/IMPI/bin/greasy $file" >> commands_to_execute_greasy.txt
done