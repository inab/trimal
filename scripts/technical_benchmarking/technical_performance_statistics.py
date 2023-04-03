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

variable_labels = {
    "repetition": "Repetitions",
    "method": "Method",
    "alignment_size": "Size (bytes)",
    "alignment_size_mb": "Size (MB)",
    "alignment_seqs": "Sequences",
    "alignment_cols": "Columns",
    "alignment_res": "Residues",
    "alignment_res_millions": "Residues (millions)",
    "trimmed_alignment_cols": "Trimmed columns",
    "user_time": "Time (s)",
    "system_time": "System time (s)",
    "percent_cpu": "CPU (%)",
    "max_resident_set_size": "Max memory (KB)",
    "max_resident_set_size_gb": "Max memory (GB)",
    "exit_status": "Exit status"
}


def generate_plot(plot_type, df, x_axis, y_axis, hue, log_x, log_y, title, show_plot, save_plot):
  if plot_type == "stripplot":
    if y_axis == "alignment_res":
      alignment_res_ordered = np.sort(df["alignment_res"].unique())[::-1]
      tech_plot = sns.stripplot(data=df, x=x_axis, y=y_axis,
                      hue=hue, orient="h", order=list(alignment_res_ordered))
      tech_plot.set_yticklabels(["{:,.0f}".format(label)
                      for label in alignment_res_ordered])
    else:
      tech_plot = sns.stripplot(data=df, x=x_axis, y=y_axis,
                        hue=hue, orient="h")
  elif plot_type == "lineplot":
    tech_plot = sns.lineplot(data=df, x=x_axis, y=y_axis, hue=hue)
  elif plot_type == "violinplot":
    tech_plot = sns.violinplot(data=df, x=x_axis, y=y_axis, hue=hue) 
  elif plot_type == "boxplot":
    tech_plot = sns.boxplot(data=df, x=x_axis, y=y_axis, hue=hue) 
  elif plot_type == "histplot":
    tech_plot = sns.histplot(data=df, x=x_axis, hue=hue, multiple="stack") 

  tech_plot.get_legend().set_title(title)
  plt.xlabel(variable_labels[x_axis])
  if y_axis: plt.ylabel(variable_labels[y_axis])
  if log_x: plt.xscale('log')
  if log_y: plt.yscale('log')
  if show_plot: plt.show()
  if save_plot: tech_plot.get_figure().savefig(f'{plot_type}_{x_axis}_{y_axis}_{hue}', dpi = 400)


def generate_time_plots():
   return


def generate_memory_plots():
   return


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--in", dest="inFile",
                        required=True, type=str, help="Input dataset")

    parser.add_argument("-p", "--plot", dest="plot", required=True,
                        type=str, choices=["lineplot", "boxplot", "stripplot", "violinplot", "histplot"], help="Set plot type")

    parser.add_argument("--all_of", dest="all_of", required=False,
                        type=str, choices=["time", "memory"], help="Generate all plots of a variable")

    parser.add_argument("-x", "--x_axis", dest="x_axis", required=False,
                        type=str, choices=["repetition", "method", "alignment_size", "alignment_seqs",
                                           "alignment_cols", "alignment_res", "alignment_res_millions",
                                               "alignment_size_mb", "trimmed_alignment_cols",
                                               "user_time", "system_time", "percent_cpu", "max_resident_set_size",
                                               "max_resident_set_size_gb", "exit_status"], help="Set x axis variable")

    parser.add_argument("-y", "--y_axis", dest="y_axis", required=False,
                        type=str, choices=["repetition", "method", "alignment_size", "alignment_seqs",
                                           "alignment_cols", "alignment_res", "alignment_res_millions",
                                               "alignment_size_mb", "trimmed_alignment_cols", "user_time",
                                               "system_time", "percent_cpu", "max_resident_set_size", "max_resident_set_size_gb"
                                               "exit_status"], help="Set y axis variable")

    parser.add_argument("--log_x", dest="log_x", default=False,
                        action="store_true", help="Set logaritmic scale in x axis")
    
    parser.add_argument("--log_y", dest="log_y", default=False,
                        action="store_true", help="Set logaritmic scale in y axis")

    parser.add_argument("--hue", dest="hue", required=True,
                        type=str, choices=["method", "alignment_seqs", "alignment_cols"], help="Set hue variable")
    
    parser.add_argument("--title", dest="title",
                        required=False, type=str, help="Set plot title")
    
    parser.add_argument("--filter_by_method", dest="method",
                        required=False, type=str, help="Set method to filter by")
    
    parser.add_argument("--filter_by_min_sequences", dest="min_sequences",
                        required=False, type=str, help="Set minimum sequences filter")
    
    parser.add_argument("--filter_by_min_columns", dest="min_columns",
                        required=False, type=str, help="Set minimum columns filter")
    
    parser.add_argument("--show", dest="show_plot", default=False,
                        action="store_true", help="Show plot interactively")

    parser.add_argument("--save", dest="save_plot", default=False,
                        action="store_true", help="Save plot")

    args = parser.parse_args()

    if not os.path.isfile(args.inFile):
        sys.exit(("ERROR: Check input alignment file '%s'") % (args.inFile))

    if not args.all_of and not args.x_axis and not args.y_axis:
        sys.exit("ERROR: all_of or axis variables should be specified")

    df = pd.read_csv(args.inFile)

    if args.method: df = df[df["method"] == args.method]
    if args.min_sequences: df = df[df["alignment_seqs"] >= args.min_sequences]
    if args.min_columns: df = df[df["alignment_cols"] >= args.min_columns]

    df["alignment_res"] = df["alignment_seqs"] * df["alignment_cols"]
    df["alignment_res_millions"] = df["alignment_res"] * 10**(-6)
    df["alignment_size_mb"] = df["alignment_size"] * \
        10**(-6)  # max_resident_set_size is in bytes
    df["max_resident_set_size_gb"] = df["max_resident_set_size"] * \
        10**(-6)  # max_resident_set_size is in KB
    df.set_index(["file", "method", "repetition"])

    generate_plot(args.plot, df, args.x_axis, args.y_axis, args.hue, args.log_x, args.log_y,
                   args.title, args.show_plot, args.save_plot)


if __name__ == "__main__":
    sys.exit(main())
