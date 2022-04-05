#!/usr/bin/python

#
#  'technical_performance_statistics.py'
#
#   Script implemented to work with trimAl to analyze technical performance
#   statistics.
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
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
  parser = argparse.ArgumentParser()

  args = parser.parse_args()

  if os.path.exists("technical_performance.csv"):
    df = pd.read_csv("technical_performance.csv", index_col = 0)
  else:
    alignment_names = pd.read_table("technical_performance_results/alignment_filenames.txt", header = None, names=["alignment"], dtype="string")
    num_sequences = pd.read_table("technical_performance_results/alignment_seqs.txt", header = None, names=["num_sequences"], dtype="int")
    num_columns = pd.read_table("technical_performance_results/alignment_cols.txt", header = None, names=["num_columns"], dtype="int")
    file_size = pd.read_table("technical_performance_results/alignment_sizes.txt", header = None, names=["size (bytes)"], dtype="int")
    max_resident_set_size = pd.read_table("technical_performance_results/max_resident_set_size.txt", header = None, names=["max_resident_set_size (kbytes)"], dtype="int")
    percent_cpu = pd.read_table("technical_performance_results/percent_cpu.txt", header = None, names=["percent_cpu"], dtype="int")
    system_time = pd.read_table("technical_performance_results/system_times.txt", header = None, names=["system_time (s)"], dtype="float")
    user_time = pd.read_table("technical_performance_results/user_times.txt", header = None, names=["user_time (s)"], dtype="float")
    exit_status = pd.read_table("technical_performance_results/exit_status.txt", header = None, names=["exit_status"], dtype="int")

    df = pd.concat([alignment_names, num_sequences, num_columns, file_size, max_resident_set_size, percent_cpu, system_time, user_time, exit_status], axis = 1)
    df["total_residues"] = df["num_sequences"] * df["num_columns"]

  if not os.path.exists("technical_performance.csv"):
    with open("technical_performance.csv", 'w') as file:
      file.write(df.to_csv())

  df["size (megabytes)"] = df["size (bytes)"] * 10**(-6)
  df["total_residues (millions)"] = df["total_residues"] * 10**(-6)
  plt.plot("user_time (s)", "total_residues (millions)", data = df, linestyle="none", marker="o")
  plt.xlabel("Time (s)")
  plt.ylabel("Residues (millions)")
  plt.title("Read and write of alignment with trimAl")
  plt.show()


if __name__ == "__main__":
  sys.exit(main())
