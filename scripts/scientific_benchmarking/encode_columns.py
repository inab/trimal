#!/usr/bin/python

#
#  'encode_columns.py'
#
#   Script implemented to encode alignment columns
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


import argparse
from ctypes import alignment
import hashlib
import numpy as np
import os
import sys

from Bio import AlignIO


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in", dest = "inFile", required = True, type = str, help = "Input alignment")
    parser.add_argument("-f", "--format", dest = "inFormat", default = "fasta", \
    type = str, choices = ["clustal", "fasta-m10", "fasta", "phylip-relaxed", \
    "phylip-sequential", "phylip", "nexus"],help = "Set input alignment format")

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

    alignment = AlignIO.read(args.inFile, format = args.inFormat)
    align_array = np.array([list(rec) for rec in alignment], order="F")
    align_array = np.transpose(align_array)
    for col in align_array:
        hash_obj = hashlib.md5(col)
        print(hash_obj.hexdigest())


if __name__ == "__main__":
  sys.exit(main())