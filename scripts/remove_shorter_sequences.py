#!/usr/bin/python

#
# 'remove_shorter_sequences.py'
#
#   Script implemented to explore future functionalities of trimAl. The script
#   analyzes the length of each sequence and remove those shorter than a given
#   length set by the user
#
#   [2015] S. Capella-Gutierrez - scapella@crg.es
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

  parser.add_argument("-m", "--min", dest = "minLen", default = 1, type = int,
    help = "Set a minimum sequence length to keep it in the output alignment")

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

  ofile = open(args.outFile, "w") if args.outFile else sys.stdout
  for record in AlignIO.read(args.inFile, format = args.inFormat):
    sequence_id = record.id if not args.keepHeader else record.description
    sequence = str(record.seq)

    length = len(sequence)
    valid = len([ps for ps in range(length) if sequence[ps] != args.gapSymbol])

    if valid >= args.minLen:
      print >> ofile, (">%s\n%s") % (sequence_id, sequence)
    elif args.verbose:
      msg =  ("INFO: Sequence '%s' has been removed. Shorter ") % (sequence_id)
      msg += ("(%d) than min. sequence length (%d)") % (valid, args.minLen)
      print >> sys.stderr, msg
      sys.stderr.flush()
  ofile.close()
