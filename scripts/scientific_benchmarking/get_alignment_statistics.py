#!/usr/bin/python

#
# 'get_alignment_statistics.py'
#
#   Script implemented to work with trimAl to analyze statistics of
#   aligned sequences.
#
#   [2022] Nicolás Díaz Roussel - nicolas.diazroussel@bsc.es
#
#   this script is free software: you can redistribute it and/or modify it under
#   the terms of the GNU General Public License as published by the Free
#   Software Foundation, the last available version.
#
#   this script is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
#   more details on <http://www.gnu.org/licenses/>
#
from Bio import AlignIO, SeqIO
from ete3 import Tree
import argparse
import sys
import os
import subprocess
import pathlib


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in", dest="inFile",
                        required=True, type=str, help="Input alignment")

    parser.add_argument("-rf", "--ref_tree", dest="ref_tree",
                        required=False, type=str, help="Input reference tree")

    parser.add_argument("-t", "--type", dest="residueType", required=True,
                        type=str, choices=["AA", "DNA"], help="Residue type")

    parser.add_argument("--taxon", dest="taxon", required=False, type=str,
                        choices=["Bacteria", "Eukaryotes", "Fungi"], help="Taxon")

    parser.add_argument("--single", dest="single", required=False, type=str,
                        choices=["all_gaps_indets_sequence"], help="Single statistics to calculate")

    args = parser.parse_args()
    msa_filename, msa_file_extension = os.path.splitext(args.inFile)

    if not os.path.isfile(args.inFile) or (msa_file_extension != ".fa" and msa_file_extension != ".fasta"):
        sys.exit(("Error loading file '%s'") % (args.inFile))

    if (msa_file_extension != ".fa" and msa_file_extension != ".fasta"):
        sys.exit(("Error loading file '%s'. Supported extensions are .fa and .fasta") % (
            args.inFile))

    msa_basename_splits =  os.path.basename(msa_filename).split('.')
    if len(msa_basename_splits) != 3:
        sys.exit(("Error loading file '%s'."
                  " Filename should follow the format problem_number.msa_tool.trimming_method.fa") % (args.inFile))

    if args.ref_tree:
        ref_tree_filename, ref_tree_file_extension = os.path.splitext(
            args.ref_tree)
        if not os.path.isfile(args.ref_tree) or ref_tree_file_extension != ".nwk":
            sys.exit(("ERROR: Check input reference tree file '%s'") %
                     (args.ref_tree))

    if args.single:
        if args.single == "all_gaps_indets_sequence":
            print(are_all_gaps_indets(args.inFile, args.residueType))
    else:
        if not args.taxon:
            sys.exit(("ERROR: Set taxon [--taxon]"))
        generate_statistics(args.inFile, args.residueType,
                            args.taxon, args.ref_tree)


def are_all_gaps_indets(file, residue_type):
    gap_symbol = "-"
    indet = "X" if residue_type == "AA" else "N"
    for seqRecord in SeqIO.parse(file, format="fasta"):
        if all((residue == gap_symbol or residue == indet) for residue in seqRecord.seq):
            return True

    return False


def generate_statistics(msa_filepath, residue_type, taxon, ref_tree_filepath):
    trimal_dir = pathlib.Path(__file__).parent.parent.parent.resolve()
    trimal_bin = f'{trimal_dir}/bin/trimal'
    msa_basename = os.path.basename(msa_filepath)
    msa_id = f'{residue_type}_{taxon}_{msa_basename}'
    problem_num = msa_basename.split('.')[0][-4:]
    msa_tool = msa_basename.split('.')[1]
    msa_filter_tool = msa_basename.split('.')[2]
    temp_trimal_sgc_filename = f'temp_trimal_sgc_{msa_id}.txt'
    with open(temp_trimal_sgc_filename, "w") as temp_trimal_sgc_output:
        subprocess.run([trimal_bin, "-in", msa_filepath, "-sgc"],
                       stdout=temp_trimal_sgc_output)

    temp_number_blocks_filename = f'temp_number_blocks_output_{msa_id}.txt'
    with open(temp_number_blocks_filename, "w") as temp_number_blocks_output:
        subprocess.run(["python3", f'{trimal_dir}/scripts/set_manual_boundaries.py', "-i",
                        temp_trimal_sgc_filename, "--inner_blocks",
                        "--total_blocks", "--block_coordinates",
                        "--min_gapscore_allowed", "1", "--max_blocks", "10000"],
                       stdout=temp_number_blocks_output)

    number_blocks = -1
    left_block_column = -1
    right_block_column = -1
    with open(temp_number_blocks_filename, "r") as temp_number_blocks_output:
        for line in temp_number_blocks_output.readlines():
            if "## Blocks" in line:
                number_blocks = line.split()[-1]
            elif "## Left column" in line:
                left_block_column = line.split()[-1]
            elif "## Right column" in line:
                right_block_column = line.split()[-1]

    msa = AlignIO.read(msa_filepath, format="fasta")
    number_columns = msa.get_alignment_length()
    number_sequences = len(msa)

    temp_trimal_sgt_filename = f'temp_trimal_sgt_{msa_id}.txt'
    with open(temp_trimal_sgt_filename, "w") as temp_trimal_sgt_output:
        subprocess.run([trimal_bin, "-in", msa_filepath, "-sgt"],
                       stdout=temp_trimal_sgt_output)
    avg_gaps = -1
    with open(temp_trimal_sgt_filename, "r") as temp_trimal_sgt_output:
        for line in temp_trimal_sgt_output.readlines():
            if "## AverageGaps" in line:
                avg_gaps = line.split()[-1]
                break

    temp_trimal_sident_filename = f'temp_trimal_sident_{msa_id}.txt'
    with open(temp_trimal_sident_filename, "w") as temp_trimal_sident_output:
        subprocess.run([trimal_bin, "-in", msa_filepath, "-sident"],
                       stdout=temp_trimal_sident_output)
    avg_seq_identity = -1
    with open(temp_trimal_sident_filename, "r") as temp_trimal_sident_output:
        for line in temp_trimal_sident_output.readlines():
            if "## AverageIdentity" in line:
                avg_seq_identity = line.split()[-1]
                break

    subprocess.run(["rm", temp_trimal_sgc_filename, temp_number_blocks_filename,
                    temp_trimal_sgt_filename, temp_trimal_sident_filename])

    rf_distance = -1
    if ref_tree_filepath:
        t1 = Tree(f'{msa_filepath}.treefile')
        t2 = Tree(ref_tree_filepath)
        try:
            rf_distance, max_parts, common_attrs, edges1, edges2, discard_t1, discard_t2 = t1.robinson_foulds(
                t2, unrooted_trees=True)
        except:
            # write error into log
            print()
    
    print(f'{residue_type},{taxon},{problem_num},{msa_tool},{msa_filter_tool},{number_columns},{number_sequences},{number_blocks},{left_block_column},{right_block_column}',
          f'{avg_gaps},{avg_seq_identity},{rf_distance}')


def generate_comparison_statistics():
    # "$columns_raw_MSA,$removed_columns,$percent_conserved_columns,$avg_gaps_diff_weighted,$avg_seq_identity_diff_weighted,$RF_distance_diff"
    print()


if __name__ == "__main__":
    sys.exit(main())


