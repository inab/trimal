#!/usr/bin/python

#
#  'technical_performance_statistics.py'
#
#   Script implemented to work with trimAl to analyze technical performance
#   statistics.
#
#   [2022] Nicolás Díaz Roussel - nicolas.diazroussel@bsc.es
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

  df_none = df.loc[df['method'] == 'None', :]
  sns.lineplot(y=df_none["alignment_size"] * 10**(-6), x=df_none["user_time"], marker="o")
  plt.xlabel("Time (s)")
  plt.ylabel("Size (megabytes)")
  plt.title("Input/output of MSA with trimAl")
  plt.xscale('log')
  plt.yscale('log')
  plt.show()

    for file in pd.unique(df["file"]):
    for repetition in pd.unique(df.loc[df["file"] == file, "repetition"]):
      df.loc[(df["file"] == file) & (df['repetition'] == repetition), "R/W_memory"] = df.loc[(df['file'] == file) & (df['method'] == "None") & (df['repetition'] == repetition), "max_resident_set_size"].values[0]
      print(file + "; " + str(repetition) + "; " + str(df.loc[(df["file"] == file) & (df['repetition'] == repetition), "R/W_memory"].values))

  '''

  #df = pd.read_csv("tech_bench_greasy_chunk1_from_5GB_to_20GB_23589053.csv")
  #df = pd.read_csv("tech_bench_greasy_chunk1_from_5GB_to_20GB_clean.csv")
  df = pd.read_csv("tech_bench_greasy_chunk1_up_to_20GB.csv")

  #df["size (megabytes)"] = df["size (bytes)"] * 10**(-6)
  #df["total_residues (millions)"] = df["total_residues"] * 10**(-6)
  #df["max_resident_set_size (gbytes)"] = df["max_resident_set_size (kbytes)"] * 10**(-6)
  print(df.head())
  print(df.describe())

  # keep the msas with more sequences instead of less sequences and more columns

  df.set_index(["file", "method", "repetition"])

  # max_resident_set_size is originally in kilobytes!!
  #df["percentage_mem_usage"] = df["max_resident_set_size"] / df["R/W_memory"]

  #df = df.loc[df["alignment_seqs"] >= df["alignment_cols"], :]
  #df = df.loc[df["file"].str.contains("example.004"), :]
  
  #df = df.loc[df["method"] == 'automated2', :]
  df["alignment_res"] = df["alignment_seqs"] * df["alignment_cols"]
  #df = df.loc[df["alignment_res"] < 5*10**8, :]
  #df["alignment_size"] = df["alignment_size"].astype("category")
  #order_sizes = np.sort(pd.unique(df["alignment_size"]))[::-1]
  #print(order_sizes)
  
  # añadir hue file para strict
  # para que no salgan disntintos tamaños con distintas columnas habría que escoger sólo una cantidad de colimnas por tamaño o de secuencias
  #ax = sns.pointplot(data=df, y="alignment_res", x="user_time", join=False, orient='h', hue="method")


  ax = sns.lineplot(data=df, y="alignment_size", x="user_time", hue="method", style="method", markers=False, dashes=False)

  ax.get_legend().set_title('')
  plt.xlabel("Time (s)")
  plt.ylabel("Size (bytes)")
  plt.title("Trimming of MSA with trimAl")
  plt.xscale('log')
  plt.yscale('log')
  plt.show()

  ax = sns.lineplot(data=df, y="alignment_seqs", x="user_time", hue="method", style="method", markers=False, dashes=False)

  ax.get_legend().set_title('')
  plt.xlabel("Time (s)")
  plt.ylabel("Sequences")
  plt.title("Trimming of MSA with trimAl")
  plt.xscale('log')
  plt.yscale('log')
  plt.show()

  ax = sns.lineplot(data=df, y="user_time", x="alignment_res", hue="method", style="method", markers=False, dashes=False)

  ax.get_legend().set_title('')
  plt.ylabel("Time (s)")
  plt.xlabel("Residues")
  plt.title("Trimming of MSA with trimAl")
  plt.xscale('log')
  plt.yscale('log')
  plt.show()
  

  ax = sns.lineplot(data=df, y="alignment_size", x="percentage_mem_usage", hue="method", style="method", markers=False, dashes=False)
  plt.xlabel("Percentage memory")
  plt.ylabel("Size (bytes)")
  plt.title("Trimming of MSA with trimAl")

  ax.get_legend().set_title('')
  plt.xscale('log')
  plt.yscale('log')
  plt.show()

  ax = sns.violinplot(data=df, y="percentage_mem_usage", x="method", hue="method")
  plt.xlabel("Method")
  plt.ylabel("Percentage memory")
  plt.title("Trimming of MSA with trimAl")
  ax.get_legend().set_title('')
  plt.show()



  
  
  return
  df.loc[:, "perc_conserved_alignment"] = df.loc[:, "trimmed_alignment_cols"] / df.loc[:, "alignment_cols"]
  df_clean = df.loc[(df["perc_conserved_alignment"] <= 1), :] # ver qué pasa en estos casos
  df_clean = df_clean.loc[df_clean["method"] != "None", :]
  #print(df.loc[df["method"] == "None", "perc_conserved_alignment"])
  ax = sns.lineplot(y=df_clean["perc_conserved_alignment"], x=df_clean["user_time"], hue=df_clean["method"], style=df_clean["method"], markers=False, dashes=False)
  plt.xlabel("Time (s)")
  plt.ylabel("% of conserved columns")
  plt.title("Trimming of MSA with trimAl")

  ax.get_legend().set_title('')
  #plt.xscale('log')
  plt.show()

  ax = sns.scatterplot(y=df_clean["perc_conserved_alignment"], x=df_clean["user_time"], hue=df_clean["method"], style=df_clean["method"])
  plt.xlabel("Time (s)")
  plt.ylabel("% of conserved columns")
  plt.title("Trimming of MSA with trimAl")

  ax.get_legend().set_title('')
  #plt.xscale('log')
  plt.show()


if __name__ == "__main__":
  sys.exit(main())
