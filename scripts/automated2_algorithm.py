#!/usr/bin/python


import argparse
import sys
import os
import subprocess


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
      str, help = "Input alignment")

  parser.add_argument("-a", "--algorithm", dest = "algorithm", required = True, type = str, choices = ["gaps", "similarity", "combined"], help = "Method")
  parser.add_argument("--min_columns", dest = "minColumns", default = 200, required = False, type = int,  help = "Minimin number of columns to keep in the alignment")
  parser.add_argument("--min_percentage_columns", dest = "minPercentageColumns", default = 0.4, required = False, type = float,  help = "Minimin number of columns to keep in the alignment")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
      sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

  # anlayse which columns to trim based on gap score until percentage conserved columns 0.4 or num columns = 200; then trim with the flag -cols

  min_columns_to_conserve = args.minColumns
  min_percentage_columns_to_conserve = args.minPercentageColumns

  if args.algorithm == "gaps":
    p1 = subprocess.Popen(["trimal/bin/trimal", "-in", args.inFile, "-sgc"], stdout = subprocess.PIPE)
    p2 = subprocess.Popen(["tail", "-n", "+8"], stdin = p1.stdout, stdout = subprocess.PIPE)
    result = subprocess.run(["awk", "{print $3}"], stdin = p2.stdout, stdout = subprocess.PIPE)
    gaps_scores_strings = result.stdout.decode("utf-8").strip().split('\n')
    gaps_scores = [float(score) for score in gaps_scores_strings]
    sorted_indeces_gap_scores = sorted(range(len(gaps_scores)), key=lambda k: gaps_scores[k])
    alignment_length = len(gaps_scores)

    columns_to_conserve = int(min_percentage_columns_to_conserve * alignment_length) \
      if min_percentage_columns_to_conserve * alignment_length >= min_columns_to_conserve \
      else min(min_columns_to_conserve, alignment_length)

    columns_to_remove = alignment_length - columns_to_conserve
    #filename = os.path.splitext(os.path.basename(args.inFile))[0] + "_trimAl_automated2.fa"
    subprocess.run(["trimal/bin/trimal", "-in", args.inFile, "-selectcols", "{",
        ','.join(str(e) for e in sorted_indeces_gap_scores[0:columns_to_remove]), "}"])
  
  elif args.algorithm == "similarity":
    p1 = subprocess.Popen(["trimal/bin/trimal", "-in", args.inFile, "-ssc"], stdout = subprocess.PIPE)
    p2 = subprocess.Popen(["tail", "-n", "+9"], stdin = p1.stdout, stdout = subprocess.PIPE)
    result = subprocess.run(["awk", "{print $2}"], stdin = p2.stdout, stdout = subprocess.PIPE)
    similarity_scores_strings = result.stdout.decode("utf-8").strip().split('\n')
    similarity_scores = [float(score) for score in similarity_scores_strings]
    sorted_indeces_similarity_scores = sorted(range(len(similarity_scores)), key=lambda k: similarity_scores[k])
    alignment_length = len(similarity_scores)

    columns_to_conserve = int(min_percentage_columns_to_conserve * alignment_length) \
      if min_percentage_columns_to_conserve * alignment_length >= min_columns_to_conserve \
      else min(min_columns_to_conserve, alignment_length)
    
    columns_to_remove = alignment_length - columns_to_conserve
    subprocess.run(["trimal/bin/trimal", "-in", args.inFile, "-selectcols", "{",
        ','.join(str(e) for e in sorted_indeces_similarity_scores[0:columns_to_remove]), "}"])

  elif args.algorithm == "combined":
    p1 = subprocess.Popen(["trimal/bin/trimal", "-in", args.inFile, "-sgc"], stdout = subprocess.PIPE)
    p2 = subprocess.Popen(["tail", "-n", "+8"], stdin = p1.stdout, stdout = subprocess.PIPE)
    result = subprocess.run(["awk", "{print $3}"], stdin = p2.stdout, stdout = subprocess.PIPE)
    gaps_scores_strings = result.stdout.decode("utf-8").strip().split('\n')
    gaps_scores = [float(score) for score in gaps_scores_strings]
    p1 = subprocess.Popen(["trimal/bin/trimal", "-in", args.inFile, "-ssc"], stdout = subprocess.PIPE)
    p2 = subprocess.Popen(["tail", "-n", "+9"], stdin = p1.stdout, stdout = subprocess.PIPE)
    result = subprocess.run(["awk", "{print $2}"], stdin = p2.stdout, stdout = subprocess.PIPE)
    similarity_scores_strings = result.stdout.decode("utf-8").strip().split('\n')
    similarity_scores = [float(score) for score in similarity_scores_strings]
    combined_scores = [gaps_scores[i] * similarity_scores[i] for i in range(len(gaps_scores))]

    sorted_indeces_combined_scores = sorted(range(len(combined_scores)), key=lambda k: combined_scores[k])
    alignment_length = len(combined_scores)

    columns_to_conserve = int(min_percentage_columns_to_conserve * alignment_length) \
      if min_percentage_columns_to_conserve * alignment_length >= min_columns_to_conserve \
      else min(min_columns_to_conserve, alignment_length)
    
    columns_to_remove = alignment_length - columns_to_conserve
    subprocess.run(["trimal/bin/trimal", "-in", args.inFile, "-selectcols", "{",
        ','.join(str(e) for e in sorted_indeces_combined_scores[0:columns_to_remove]), "}"])


if __name__ == "__main__":
  sys.exit(main())
