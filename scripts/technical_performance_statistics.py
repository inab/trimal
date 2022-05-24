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

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def main():

  '''
  alignment_names = pd.read_table("technical_performance_results/alignment_filenames.txt", header = None, names=["alignment"], dtype="string")
  num_sequences = pd.read_table("technical_performance_results/alignment_seqs.txt", header = None, names=["num_sequences"], dtype="int")
  num_columns = pd.read_table("technical_performance_results/alignment_cols.txt", header = None, names=["num_columns"], dtype="int")
  file_size = pd.read_table("technical_performance_results/alignment_sizes.txt", header = None, names=["size (bytes)"], dtype="int")
  max_resident_set_size = pd.read_table("technical_performance_results/max_resident_set_size.txt", header = None, names=["max_resident_set_size (kbytes)"], dtype="int")
  percent_cpu = pd.read_table("technical_performance_results/percent_cpu.txt", header = None, names=["percent_cpu"], dtype="int")
  system_time = pd.read_table("technical_performance_results/system_times.txt", header = None, names=["system_time (s)"], dtype="float")
  user_time = pd.read_table("technical_performance_results/user_times.txt", header = None, names=["user_time (s)"], dtype="float")
  exit_status = pd.read_table("technical_performance_results/exit_status.txt", header = None, names=["exit_status"], dtype="int")
  '''

  df = pd.read_csv("tech_bench_results.csv")

  #df["size (megabytes)"] = df["size (bytes)"] * 10**(-6)
  #df["total_residues (millions)"] = df["total_residues"] * 10**(-6)
  #df["max_resident_set_size (gbytes)"] = df["max_resident_set_size (kbytes)"] * 10**(-6)
  df_none = df.loc[df['method'] == 'None', :]
  sns.lineplot(x=df_none["alignment_size"] * 10**(-6), y=df_none["user_time"], marker="o")
  #plt.plot(df["user_time"], df["alignment_size"] * 10**(-6) , linestyle="none", marker="o", hue = df["method"])
  plt.ylabel("Time (s)")
  plt.xlabel("Size (megabytes)")
  plt.title("Input/output of MSA with trimAl")
  plt.show()

  df_methods = df.loc[df['method'] != 'None', :]
  ax = sns.lineplot(x=df_methods["alignment_size"] * 10**(-6), y=df_methods["user_time"], hue=df_methods["method"], marker="o")
  plt.ylabel("Time (s)")
  plt.xlabel("Size (megabytes)")
  plt.title("Trimming of MSA with trimAl")
  ax.get_legend().set_title('')
  plt.show()


if __name__ == "__main__":
  sys.exit(main())
