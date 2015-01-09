#!/usr/bin/python
from Bio import SeqIO
from string import upper
import numpy as np
import argparse
import sys
import os

codon_table = {
  ## Leucine (Leu)
  "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
  "CTN": "L",
  ## Serine (Ser)
  "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
  "TCN": "S",
  ## Arginine (Arg)
  "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
  "CGN": "R",
  ## Proline (Pro)
  "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
  "CCN": "P",
  ## Glycine (Gly)
  "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
  "GGN": "G",
  ## Alanine (Ala)
  "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
  "GCN": "A",
  ## Valine (Val)
  "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
  "GTN": "V",
  ## Threonine (Thr)
  "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
  "ACN": "T",
  ## Isoleucine (Ile)
  "ATT": "I", "ATC": "I", "ATA": "I",
  ## Phenylalanine (Phe)
  "TTT": "F", "TTC": "F",
  ## Tyrosine (Tyr)
  "TAT": "Y", "TAC": "Y",
  ## Cysteine (Cys)
  "TGT": "C", "TGC": "C",
  ## Histidine (His)
  "CAT": "H", "CAC": "H",
  ## Glutamine (Gln)
  "CAA": "Q", "CAG": "Q",
  ## Aspartic acid (Asp)
  "GAT": "D", "GAC": "D",
  ## Glutamic acid (Glu)
  "GAA": "E", "GAG": "E",
  ## Lysine (Lys)
  "AAA": "K", "AAG": "K",
  ## Asparagine (Asn)
  "AAT": "N", "AAC": "N",
  ## Tryptophan (Trp)
  "TGG": "W", 
  ## Methionine (M), Start
  "ATG": "M",

  ## Stop codons
  "TGA": "U", ## Selenocysteine (Sel)
  "TAG": "O", ## Pyrrolysine (Pyl)
  "TAA": "X",

  ## Additional characters
  "NNN": "X",
}

stop_codons = {
  "TGA": "U", ## Selenocysteine (Sel)
  "TAG": "O", ## Pyrrolysine (Pyl)
  "TAA": "X",
}

def _split(seq, length = 80):
  return "\n".join([seq[i:i + length] for i in range(0, len(seq), length)])

if __name__ == "__main__":

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input Codon alignment")

  parser.add_argument("-o", "--out", dest = "outFile", default = None, type = \
    str, help = "Set output file")

  parser.add_argument("-l", "--log", dest = "logFile", default = None, type = \
    str, help = "Set output log file")

  parser.add_argument("-w", "--windows_size", dest = "wSize", default = 2, \
    type = int, help = "Set how many columns should be analyzed before/after "
    + "each position")

  parser.add_argument("-f", "--format", dest = "inFormat", default = "fasta", \
    type = str, choices = ["clustal", "fasta-m10", "fasta", "phylip-relaxed", \
    "phylip-sequential", "phylip", "nexus"],help = "Set input alignment format")

  parser.add_argument("-g", "--gap_symbol", dest = "gapSymbol", default = '-', \
    type = str, help = "Define the gap symbol used in the input alignment")

  parser.add_argument("--discard_gaps", dest = "noGaps", action = "store_true",
    default = False, help = "Discard any column containing gaps prior to any "
    + "analysis")

  parser.add_argument("-v", "--verbose", dest = "verbose", default = True,
    action = "store_false", help = "Deactivate verbosity")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input CDS alignment file '%s'") % (args.inFile))

  if args.verbose:
    logFile = open(args.logFile, "w") if args.logFile else sys.stderr

  if args.wSize < 0:
    sys.exit(("ERROR: Check windows size '%s'") % (str(args.winSize)))

  codon_table.setdefault(("%s") % (args.gapSymbol * 3), args.gapSymbol)

  alignment_nt, alignment, incongruences = {}, {}, {}
  order, aligLength = [], 0

  incongruentCodonsCols = set()  
  for record in SeqIO.parse(args.inFile, args.inFormat):
    seq = str(record.seq)
    if record.id in alignment:
      sys.exit(("ERROR: Duplicated entry '%s'") % (record.id))
    if aligLength and aligLength != len(seq):
      sys.exit("ERROR: Check alignment - it seems sequences are unaligned")

    aligLength = len(seq)
    order.append(record.id)
    alignment_nt.setdefault(record.id, seq)

    translated = ""
    for pos in range(3, aligLength+1, 3):
      codon = seq[pos-3:pos]
      if not codon in codon_table:
        incongruences.setdefault(record.id, set()).add((pos/3)-1)
        incongruentCodonsCols.add((pos/3)-1)
      translated += codon_table[codon] if codon in codon_table else "X"
    alignment.setdefault(record.id, translated)
   
  ## Store the initial alignment size
  initialLength = aligLength/3
  incongruentCodons = len(incongruentCodonsCols)

  selected = sorted(set(range(initialLength)) - incongruentCodonsCols)
  aligLength = len(selected)

  if args.verbose and incongruences:
    print >> logFile, ("## Incongruences in Codons > AAs translation")
    for record in sorted([(len(incongruences[s]), s) for s in incongruences]):
      print >> logFile, ("%-12s\t%8d") % (record[1], record[0])

    msg = "Initial Alignment Size:"
    print >> logFile, ("\n## Stats\n%-58s\t%8d") % (msg, initialLength * 3)
    msg = "Columns containing at least one untranslatable codons:"
    ratio = incongruentCodons/float(initialLength)
    print >> logFile, ("%-58s\t%8d\t%.4f") % (msg, incongruentCodons * 3, ratio)
    msg = "Alignment Size after trimming:"
    ratio = aligLength/float(initialLength)
    print >> logFile, ("%-58s\t%8d\t%.4f") % (msg, aligLength * 3, ratio)
 
  ## If requested, remove any column containing gaps prior to any analysis
  discardedGappyCols = set()
  if args.noGaps:
    for pos in selected:
      if args.gapSymbol in [alignment[seq][pos] for seq in alignment]:
        discardedGappyCols.add(pos)
    discardedGappy = len(discardedGappyCols)

  ## Update which columns could be used in downstream analyses
  toRemove = discardedGappyCols | incongruentCodonsCols
  selected = sorted(set(range(initialLength)) - toRemove)
  aligLength = len(selected)

  nonConservedNeighboursCols = set()
  
  ## Analyze alignment extremes: right
  size = args.wSize * 2
  upper_end = size + 1
  for pos in range(args.wSize):
    ## Check if the given column is conserved or not
    #~ if len(set([alignment[seq][selected[pos]] for seq in alignment])) == 1:
      #~ continue

    ## Check surrounding columns to see whether all of them are conserved or not
    conserve = True
    for col in range(pos) + range(pos+1, upper_end):
      if len(set([alignment[seq][selected[col]] for seq in alignment])) != 1:
        conserve = False
        break

    if not conserve:
      nonConservedNeighboursCols.add(selected[pos])

  ## Analyze alignment extremes: left
  lower_start = aligLength - size - 1
  for pos in range(aligLength - args.wSize, aligLength):
    ## Check if the given column is conserved or not
    #~ if len(set([alignment[seq][selected[pos]] for seq in alignment])) == 1:
      #~ continue

    ## Check surrounding columns to see whether all of them are conserved or not
    conserve = True
    for col in range(lower_start, pos) + range(pos+1, aligLength):
      if len(set([alignment[seq][selected[col]] for seq in alignment])) != 1:
        conserve = False
        break

    if not conserve:
      nonConservedNeighboursCols.add(selected[pos])

  ## Analyze the rest of the alignment
  for pos in range(args.wSize, aligLength - args.wSize):

    ## Check whether the current column is fully conserved or not -
    ## In case of being fully conserved, move to next column
    #~ if len(set([alignment[seq][selected[pos]] for seq in alignment])) == 1:
      #~ continue

    ## Check surrounding columns to see whether all of them are conserved or not
    conserve = True
    for col in range(pos - args.wSize, pos) + range(pos+1, pos+1 + args.wSize):
      if len(set([alignment[seq][selected[col]] for seq in alignment])) != 1:
        conserve = False
        break

    if not conserve:
      nonConservedNeighboursCols.add(selected[pos])
  nonConservedNeighbours = len(nonConservedNeighboursCols)

  ## Update with which columns should be removed and which ones kept
  toRemove |= nonConservedNeighboursCols
  selected = sorted(set(range(initialLength)) - toRemove)

  ## Print some report
  if args.noGaps and args.verbose:
    if not incongruences:
      msg = "Initial Alignment Size:"
      print >> logFile, ("## Stats\n%-58s\t%8d") % (msg, initialLength * 3)
    msg = "Columns containing at least 1 gaps:"
    ratio = discardedGappy/float(initialLength)
    print >> logFile, ("%-58s\t%8d\t%.4f") % (msg, discardedGappy * 3, ratio)
    msg = "Alignment Size after trimming:"
    ratio = aligLength/float(initialLength)
    print >> logFile, ("%-58s\t%8d\t%.4f") % (msg, aligLength * 3, ratio)
    
  if args.verbose:
    if not discardedGappyCols and not incongruences:
      msg = "\nInitial Alignment Size:"
      print >> logFile, ("## Stats\n%-58s\t%8d") % (msg, initialLength * 3)
      
    msg = "Columns with non-conserved neighbours:"
    r = nonConservedNeighbours/float(aligLength)
    print >> logFile,("%-58s\t%8d\t%.4f") % (msg, nonConservedNeighbours * 3, r)
    msg = "Alignment Size after trimming:"
    final = (aligLength-nonConservedNeighbours)
    ratio = final/float(initialLength)
    print >> logFile, ("%-58s\t%8d\t%.4f") % (msg, final * 3, ratio)
   
    if args.logFile:
      output = ",".join(map(str, sorted(toRemove)))
      print >> logFile, ("## Discarded Columns\t%s") % (output)

  ofile = open(args.outFile, "w") if args.outFile else sys.stdout
  if final > 0:
    for seqName in order:
      output = "".join([alignment_nt[seqName][3*pos:3*(pos+1)] for pos in selected])
      print >> ofile, (">%s\n%s") % (seqName, _split(output))
  ofile.close()
