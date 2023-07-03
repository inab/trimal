#!/usr/bin/python3

#
# 'get_sequneces_gaps.py'
#
#   Script implemented to obtain the sequences index for those seuqneces 
#   exceding a minimum gaps' ratio threshold.
#
#   [2020] S. Capella-Gutierrez - salvador.capella@bsc.es
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
from Bio import AlignIO
import argparse
import sys
import os

if __name__ == "__main__":

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input alignment")

  parser.add_argument("-o", "--out", dest = "outFile", default = None, type = \
    str, help = "Set output file. It will be generated into FASTA format")

  parser.add_argument("-f", "--format", dest = "inFormat", default = "fasta", \
    type = str, choices = ["clustal", "fasta-m10", "fasta", "phylip-relaxed", \
    "phylip-sequential", "phylip", "nexus"],help = "Set input alignment format")

  parser.add_argument("-g", "--gap_symbol", dest = "gapSymbol", default = '-', \
    type = str, help = "Define the gap symbol used in the input alignment")

  parser.add_argument("--show_only_index", dest = "showIndexes", default = False, \
    action = "store_true", help = "Show only the indexes of sequences with a "
    + "gaps' ratio equal or higher than the established threshold")
    
  parser.add_argument("--threshold", dest = "gapsThreshold", default = 0.0, \
    type = float, help = "Identify sequences with a minimum of gaps' ratio")

  parser.add_argument("--keep_header", dest = "keepHeader", default = False,
    action = "store_true", help = "Keep original alignment sequence IDs indepen"
    + "dently of blank spaces on it")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

  index = 0
  indexes = []
  ofile = open(args.outFile, "w") if args.outFile else sys.stdout
  for record in AlignIO.read(args.inFile, format = args.inFormat):
    sequence_id = record.id if not args.keepHeader else record.description
    sequence = str(record.seq)

    length = len(sequence)
    valid = len([ps for ps in range(length) if sequence[ps] != args.gapSymbol])
    gaps_ratio = 1 - (valid/length)

    if gaps_ratio >= args.gapsThreshold:
      if not args.showIndexes:
        print(f'{index:4d}\t{sequence_id:30}\t{gaps_ratio:.4f}', file = ofile)
      indexes.append(index)
    index += 1

  if args.showIndexes:
    print (','.join(map(str, indexes)), file = ofile)

  ofile.close()
