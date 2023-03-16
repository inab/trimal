#!/usr/bin/python

#
# 'get_alignment_statistics.py'
#
#   Script implemented to work with trimAl to analyze statistics of
#   aligned sequences.
#
#   [2022] Nicolás Díaz Roussel - ndiazrou@bsc.es
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
from Bio import SeqIO
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
        str, help = "Input alignment")

    parser.add_argument("-t", "--type", dest = "residueType", required = True, type = \
        str, choices = ["AA", "DNA"], help = "Residue type")

    parser.add_argument("--all_gaps_indets", dest = "allGapsIndets", required = False, \
        action = "store_true", default = False, help = "Check if alignment has any "
        "sequence with only gaps and/or indets")

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))
    
    if args.allGapsIndets:
        print(are_all_gaps_indets(args.inFile, args.residueType))


def are_all_gaps_indets(file, residue_type):
    gap_symbol = "-"
    indet = "X" if residue_type == "AA" else "N"
    for seqRecord in SeqIO.parse(file, format = "fasta"):
        if all((residue == gap_symbol or residue == indet) for residue in seqRecord.seq):
            return True
    
    return False
    

if __name__ == "__main__":
  sys.exit(main())