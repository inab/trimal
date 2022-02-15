#!/usr/bin/python

#
# 'get_pre_alignment_statistics.py'
#
#   Script implemented to work with trimAl to analyze statistics of
#   pre-aligned sequences.
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

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

    min_columns = float("inf")
    max_columns = 0
    for seq in SeqIO.parse(args.inFile, format = "fasta"):
        min_columns = min(min_columns, len(seq))
        max_columns = max(max_columns, len(seq))
    
    print("## Minimum columns " + str(min_columns))
    print("## Maximum columns " + str(max_columns))


if __name__ == "__main__":
  sys.exit(main())