#!/usr/bin/python

#
# 'filter_homfam_alignments.py'
#
#   Script implemented to select the reference sequences
#   of an alignment.
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
from Bio import SeqIO
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
        str, help = "Input alignment")

    parser.add_argument("-r", "--ref", dest = "refSeqs", required = True, type = \
        str, help = "Reference sequences")

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))
    
    if not os.path.isfile(args.refSeqs):
        sys.exit(("ERROR: Check reference alignment file '%s'") % (args.refSeqs))
    
    filtered_seqs = []
    indexed_alignment = SeqIO.index(args.inFile, format = "fasta")
    
    for seq in SeqIO.parse(args.refSeqs, format = "fasta"):
        filtered_seqs.append(indexed_alignment.get(seq.id))

    filtered_seqs_filename = args.inFile + "_filtered_refs"
    SeqIO.write(filtered_seqs, filtered_seqs_filename, "fasta")



if __name__ == "__main__":
  sys.exit(main())