#!/usr/bin/python

#
# 'check_codon_alignments.py'
#
#   Script implemented to analyze resulting back-translated alignments by trimAl
#   Main idea here is to remove those codon-columns composed by only 'N'/'n' -
#   which are the symbol to indicate indeterminate nucleotides.
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

def splitSequence(seq, length = 80):
  ''' Split a given sequence contained in one line into lines of size "length"
  '''
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

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
    type = str, help = "Set the gap symbol used in the input alignment")

  parser.add_argument("--indeter_symbol", dest = "indeterSymbol", default = 'N',
    type = str, help = "Set the indetermination symbol used in the alignment")

  parser.add_argument("--keep_header", dest = "keepHeader", default = False,
    action = "store_true", help = "Keep original alignment sequence IDs indepen"
    + "dently of blank spaces on it")

  parser.add_argument("--complementary", dest = "complement", default = False,
    action = "store_true", help = "Get the complementary output alignment")

  parser.add_argument("-v", "--verbose", dest = "verbose", default = False,
    action = "store_true", help = "Activate verbosity")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

  alignment, alignment_length = {}, 0
  for record in AlignIO.read(args.inFile, format = args.inFormat):
    sequence_id = record.id if not args.keepHeader else record.description
    alignment.setdefault(sequence_id, str(record.seq))

    ## Check all sequences have the same length
    if alignment_length == 0:
      alignment_length = len(str(record.seq))
    if alignment_length != len(str(record.seq)):
      sys.exit("ERROR: Check input alignment. Sequences with different lengths")

  ## Check input alignment is multiple of 3
  if (alignment_length % 3) != 0:
    sys.exit("ERROR: Check input alignment. Its length is not multiple of 3")

  indetermination_cols = []
  indet = set([args.indeterSymbol.upper()])
  for pos in range(0, alignment_length, 3):

    onlyIndeter = True
    for col in range(pos, pos+3):
      column = set([alignment[seq][col].upper() for seq in alignment \
        if alignment[seq][col] != args.gapSymbol])
      if column ^ indet != set():
        onlyIndeter = False

    if onlyIndeter and not args.complement:
      indetermination_cols.append(pos)
    elif not onlyIndeter and args.complement:
      indetermination_cols.append(pos)

  if args.verbose and indetermination_cols:
    output = ",".join(map(str, sorted(indetermination_cols)))
    print >> sys.stderr, ("%s\t%s") % (args.inFile, output)

  ofile = open(args.outFile, "w") if args.outFile else sys.stdout
  for seq_id in alignment:
    output = "".join([alignment[seq_id][pos:pos+3] for pos in \
      range(0, alignment_length, 3) if not pos in indetermination_cols])
    print >> ofile, (">%s\n%s") % (seq_id, splitSequence(output))
  ofile.close()
