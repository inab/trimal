#!/usr/bin/python
import os
import Bio
import sys
import random
import argparse
import numpy as np
from Bio import SeqIO

def splitSequence(seq, length = 80):
  ''' Split a given sequence contained in one line into lines of size "length"
  '''
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

if __name__ == "__main__":

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input Codon alignment")

  parser.add_argument("-o", "--out", dest = "outFile", default = None, type = \
    str, help = "Set output file")

  parser.add_argument("-s", "--numb_sequences", dest = "numb_sequences", \
    default = 2, type = int, help = "Set how many sequences the output "
    + "alignment should contain")

  parser.add_argument("-r", "--numb_residues", dest = "numb_residues", \
    default = 100, type = int, help = "Set how many residues the output "
    + "alignment should contain")

  parser.add_argument("-f", "--input_format", dest = "inFormat", type = str, \
    default = "fasta", help = "Set input alignment format")

  parser.add_argument("-g", "--gap_symbol", dest = "gapSymbol", default = '-', \
    type = str, help = "Define the gap symbol used in the input/output "
    + "alignments")

  parser.add_argument("-m", "--max_attempts", dest = "attempts", default = 10, \
    type = int, help = "Define a maximum numnber of attempts when generating "
    + "a random alignment before giving it up")

  args = parser.parse_args()

  ## Check input parameters
  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

  if args.numb_sequences < 2:
    sys.exit(("ERROR: Check input sequences '%s'") % (str(args.numb_sequences)))

  if args.numb_residues < 2:
    sys.exit(("ERROR: Check input residues '%s'") % (str(args.numb_residues)))

  if args.attempts < 1:
    sys.exit(("ERROR: Check max. number of attempts '%s'") % (str(args.attempts)))

  ## Read input alignment and get some basic information from it e.g.
  ## sequences names, residues number, etc.
  algLen = -1
  alignment = {}
  for record in SeqIO.parse(args.inFile, args.inFormat):
    seq = str(record.seq)
    alignment.setdefault(record.id, seq)
    if algLen == -1:
      algLen = len(seq)
    if len(seq) != algLen:
      print("Detected Inconsistencies at Sequence's length", file = sys.stderr)

  sequences = list(alignment.keys())
  columns = list(range(algLen))

  ## Select randomly sequences and columns from the input alignment to populate
  ## the output alignment controlling there are not sequences nor columns
  ## composed only by gaps.
  
  ## This is an iterative process
  selected_seqs = []
  discarded_seqs = set()
  selected_cols = []
  discarded_cols = set()

  ## Set a counter to control how many attempts are done for generating the
  ## random alignment
  max_attempts = 0
  while True:

    while len(selected_seqs) < args.numb_sequences:
      selected = random.choice(sequences)
      if not selected in discarded_seqs:
        selected_seqs.append(selected)

    while len(selected_cols) < args.numb_residues:
      selected = random.choice(columns)
      if not selected in discarded_cols:
        selected_cols.append(selected)

    generated = {}
    for seq in selected_seqs:
      if seq in generated:
        continue
      ## We check generated sequences are not composed only by gaps
      sequence = [alignment[seq][pos] for pos in selected_cols]
      if set(sequence) - set([args.gapSymbol]) == set([]):
        discarded_seqs.add(seq)
        continue
      generated.setdefault(seq, splitSequence("".join(sequence)))

    ## We have to check there are not columns composed only by gaps
    for column in range(len(selected_cols)):
      individual_column = [generated[seq][column] for seq in generated]
      if set(individual_column) - set([args.gapSymbol]) == set([]):
        discarded_cols.add(selected_cols[column])
      
    ## We check which sequences/residues remain after controlling by those
    ## composed only by gaps
    selected_seqs = [s for s in selected_seqs if not s in discarded_seqs]
    selected_cols = [c for c in selected_cols if not c in discarded_cols]

    if len(selected_seqs) == args.numb_sequences and \
      len(selected_cols) == args.numb_residues:
      break

    max_attempts += 1
    if max_attempts == args.attempts:
      sys.exit(("ERROR: Impossible to generate random alignment after '%s' "
        + "attempts. Check configuration") % (args.attempts))
  
  ## Produce the output aligment.
  n = 1
  ofile = open(args.outFile, "w") if args.outFile else sys.stdout

  ## How to properly name output sequences including a padding to have
  ## homogeneuous ids
  padding = int(np.ceil(np.log10(args.numb_sequences)))
  if args.numb_sequences % 10 == 0:
    padding += 1
 
  for seq in selected_seqs:
    print(">seq_%s\n%s" % (str(n).zfill(padding), generated[seq]), file = ofile)
    n += 1
  ofile.close()
