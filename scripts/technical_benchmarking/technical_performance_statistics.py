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

import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def lineplot(df, show_plot):
  ax = sns.lineplot(data=df, y="alignment_size_mb", x="user_time", hue="method", style="method", markers=False, dashes=False)

  ax.get_legend().set_title('')
  plt.xlabel("Time (s)")
  plt.ylabel("Size (MB)")
  plt.title("Trimming of MSA with trimAl")
  plt.xscale('log')
  plt.yscale('log')
  if show_plot:
    plt.show()
  return

def main():

  parser = argparse.ArgumentParser()

  parser.add_argument("-i", "--in", dest = "inFile", required = True, type = \
    str, help = "Input dataset")
  
  parser.add_argument("-p", "--plot", dest = "plot", required =  True, \
    type = str, choices = ["lineplot", "boxplot", "stripplot"],help = "Set plot type")
  
  parser.add_argument("-t", "--type", dest = "type", required =  False, \
    type = str, choices = ["time", "memory"],help = "Plot all plots of a variable")
  
  parser.add_argument("-x", "--x_axis", dest = "x_axis", required =  False, \
    type = str,help = "Set x axis variable")
  
  parser.add_argument("-y", "--y_axis", dest = "y_axis", required =  False, \
    type = str,help = "Set y axis variable")
  
  parser.add_argument("--show", dest = "show_plot", default = False,
    action = "store_true", help = "Show plot interactively")
  
  parser.add_argument("--save", dest = "save_plot", default = False,
    action = "store_true", help = "Save plot")

  args = parser.parse_args()

  if not os.path.isfile(args.inFile):
    sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

  if not args.type and (not args.x_axis or not args.y_axis):
    sys.exit()
    
  df = pd.read_csv(args.inFile)
  df["alignment_res"] = df["alignment_seqs"] * df["alignment_cols"]
  df.set_index(["file", "method", "repetition"])
  
  df["alignment_res_millions"] = df["alignment_res"] * 10**(-6)
  df["alignment_size_mb"] = df["alignment_size"] * 10**(-6) # max_resident_set_size is in bytes
  df["max_resident_set_size_gb"] = df["max_resident_set_size"] * 10**(-6) # max_resident_set_size is in KB


  if args.plot == "lineplot":
    lineplot(df, args.show_plot)

  if args.plot == "stripplot":

    df = df[df["method"]=="strictplus"]
    alignment_res_ordered = np.sort(df["alignment_res"].unique())[::-1]
    print(alignment_res_ordered.shape)

    
    # sequence effect on time with respect to the same amount of residues
    ax = sns.stripplot(data=df, y="alignment_res", x="user_time", hue="alignment_seqs", orient="h", order=list(alignment_res_ordered))

    ax.get_legend().set_title('Sequences')
    plt.xlabel("Time (s)")
    plt.ylabel("Residues")
    #plt.title("Trimming of MSA with trimAl")
    plt.xscale('log')
    #plt.yscale('log')
    #for plots in ax.axes.flatten():
      #plots.yaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
      #plots.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
      #plots.yaxis.get_major_formatter().set_scientific(True)
    #plt.ticklabel_format(style='plain', axis='y')
    #ax.yaxis.set_major_formatter(mtick.StrMethodFormatter('{x:.0e}'))
    #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    print(alignment_res_ordered)
    ax.set_yticklabels(["{:,.0f}".format(label) for label in alignment_res_ordered])
    print(ax.yaxis)
    plt.show()

  return

  ax = sns.stripplot(data=df, y="alignment_res", x="user_time", hue="alignment_cols", orient="h", order=list(alignment_res_ordered))

  ax.get_legend().set_title('Columns')
  plt.xlabel("Time (s)")
  plt.ylabel("Residues")
  plt.xscale('log')
  print(alignment_res_ordered)
  ax.set_yticklabels(["{:,.0f}".format(label) for label in alignment_res_ordered])
  print(ax.yaxis)
  plt.show()


  ax = sns.boxplot(data=df, y="alignment_cols", x="user_time", hue="method", orient="h")

  ax.get_legend().set_title('')
  plt.xlabel("Time (s)")
  plt.ylabel("Columns")
  plt.title("Trimming of MSA with trimAl")
  plt.xscale('log')
  #plt.yscale('log')
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
