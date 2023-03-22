#!/usr/bin/python

#
# 'plot_alignment_statistics.py'
#
#   Script implemented to work with trimAl to plot statistics of
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
import matplotlib.pyplot as plt
import argparse
import sys
import os
import pandas as pd

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
        str, help = "Input alignment")

    parser.add_argument("--all_gaps_indets", dest = "allGapsIndets", required = False, \
        action = "store_true", default = False, help = "Plot statistics of "
        "sequence with only gaps and/or indets")

    parser.add_argument("-f", "--format", dest = "inFormat", default = "fasta", \
    type = str, choices = ["clustal", "fasta-m10", "fasta", "phylip-relaxed", \
    "phylip-sequential", "phylip", "nexus"],help = "Set input alignment format")

    parser.add_argument("--level", dest = "level", type = str, choices = ["residue",\
        "taxon", "MSA_tool", "MSA_filter_tool"], help = "Select level to plot")

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))
    
    if args.allGapsIndets:
        plot_all_gaps_indets(args.level)


def plot_all_gaps_indets(level):
    residue_type = pd.read_table("empty_seqs_statistics/residue_type.txt", header = None, names=["residue_type"], dtype="string")
    taxon = pd.read_table("empty_seqs_statistics/taxon.txt", header = None, names=["taxon"], dtype="string")
    problem_num = pd.read_table("empty_seqs_statistics/problem_num.txt", header = None, names=["problem_num"], dtype="int")
    error = pd.read_table("empty_seqs_statistics/error_problem.txt", header = None, names=["error"], dtype="string")
    msa_tools = pd.read_table("empty_seqs_statistics/MSA_tools.txt", header = None, names=["msa_tools"], dtype="string")
    msa_filter_tools = pd.read_table("empty_seqs_statistics/MSA_filter_tools.txt", header = None, names=["msa_filter_tools"], dtype="string")
    df = pd.concat([residue_type, taxon, problem_num, error, msa_tools, msa_filter_tools], axis = 1)
    if not os.path.exists("table_empty_seqs.csv"):
        with open('table_empty_seqs.csv', 'w') as file:
            file.write(df.to_csv())

    if level == "residue":
        df.groupby("residue_type").size().plot(kind='bar')
    elif level == "taxon":
        df.groupby("taxon").size().plot(kind='bar')
    elif level == "MSA_tool":
        df.groupby("msa_tools").size().plot(kind='bar')
    elif level == "MSA_filter_tool":
        df.groupby("msa_filter_tools").size().plot(kind='bar')
        
    if level != None:
        plt.savefig("empty_seqs_%s" % level, bbox_inches='tight')

if __name__ == "__main__":
  sys.exit(main())