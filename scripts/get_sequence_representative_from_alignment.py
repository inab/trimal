#!/usr/bin/python

#
# 'get_sequence_representative_from_alignment.py'
#
#   Script implemented to work with trimAl to analyze gaps statistics and decide
#   which are the boundaries in a given alignment - columns inbetween these
#   boundaries will not be removed independently of the trimming strategy
#   selected.
#
#   [2014] S. Capella-Gutierrez - scapella@crg.es
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
import numpy as np
import argparse
import sys
import os

if __name__ == "__main__":

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input alignment")

  parser.add_argument("-o", "--out", dest = "outFile", default = None, type = \
    str, help = "Set output file")

  parser.add_argument("-f", "--format", dest = "inFormat", default = "fasta", \
    type = str, choices = ["clustal", "fasta-m10", "fasta", "phylip-relaxed", \
    "phylip-sequential", "phylip", "nexus"],help = "Set input alignment format")

  parser.add_argument("-g", "--gap_symbol", dest = "gapSymbol", default = '-', \
    type = str, help = "Define the gap symbol used in the input alignment")

  parser.add_argument("--keep_header", dest = "keepHeader", default = False,
    action = "store_true", help = "Keep original alignment sequence IDs indepen"
    + "dently of blank spaces on it")

  parser.add_argument("-v", "--verbose", dest = "verbose", default = False,
    action = "store_true", help = "Activate verbosity")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

  identities, sequences = {}, {}
  for record in AlignIO.read(args.inFile, format = args.inFormat):
    current_seq = str(record.seq)
    sequence_length = len(current_seq)
    sequence_id = record.id if not args.keepHeader else record.description

    for seq in sequences:
      ## Identity score is computed considering all positions for which at least
      ## one of the sequences has a non-gap symbol
      valid_pos = [ pos for pos in range(sequence_length) if current_seq[pos] \
        != args.gapSymbol or sequences[seq][0][pos] == args.gapSymbol ]

      identical = [ pos for pos in valid_pos if sequences[seq][0][pos] == \
        current_seq[pos]]

      ratio = float(len(identical))/len(valid_pos)
      identities.setdefault(sequence_id, {}).setdefault(seq, ratio)
      identities.setdefault(seq, {}).setdefault(sequence_id, ratio)

    ## Save current sequence and  move on to the nex one
    ungapped = current_seq.replace(args.gapSymbol, "")
    sequences.setdefault(sequence_id, [current_seq, ungapped, len(ungapped)])

  selection, maxIdentity = set(), 0
  for refer in sequences:
    avg = np.average([identities[refer][seq] for seq in identities[refer]])
    if args.verbose:
      print >> sys.stderr, ("%-20s\t%.6f") % (refer, avg)
    ## Save current sequence if it has a greater identity score 
    if avg > maxIdentity:
      maxIdentity = avg
      selection = set([(sequences[refer][1], refer)])
    elif avg == maxIdentity:
      selection |= set([(sequences[refer][1], refer)])

  representative = sorted(selection, reverse = True)[0][1]
  ofile = open(args.outFile, "w") if args.outFile else sys.stdout
  print >> ofile, (">%s\n%s") % (representative, sequences[representative][1])
  ofile.close()
