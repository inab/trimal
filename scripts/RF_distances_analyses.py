#!/usr/bin/python

#
#  'parse_RF_distances.py'
#
#   Script implemented to parse RF distances data
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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("-t", "--type", dest = "residueType", required = False, type = str, choices = ["AA", "DNA"], help = "Residue type")
  parser.add_argument("--taxon", dest = "taxon", required = False, type = str, choices = ["Bacteria", "Eukaryotes", "Fungi"], help = "Taxon")
  parser.add_argument("--msa_tool", dest = "msaTool", required = False, type = str, choices = ["ClustalW", "ClustalW2", "Mafft", "Prank", "T-Coffee"], help = "MSA tool")
  parser.add_argument("--msa_filter_tool", dest = "msaFilterTool", required = False, type = str, choices = ["trimAl",], help = "MSA filter tool")


  args = parser.parse_args()

  df = pd.read_csv("RF_distances_arrays_clean.csv", index_col = 0)
  df = df.loc[df["RF_distance"] != -1, :]
  #df['RF_distance_diff'] = df['RF_distance_diff'].astype(int)
  df.loc[df["RF_distance_diff"] < 0 , 'RF_change'] = "worse"
  df.loc[df["RF_distance_diff"] > 0 , 'RF_change'] = "better"
  df.loc[df["RF_distance_diff"] == 0 , 'RF_change'] = "unchanged"

  MSA_filter_tools = pd.unique(df['MSA_filter_tool'])
  MSA_filter_tools_best = []

  for MSA_filter_tool in MSA_filter_tools:
    #RF_diff_mean = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].mean()
    #RF_diff_median = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].median()
    #RF_diff_std = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].std()
    #print("Mean, median and std RF_diff for %s are %.2f, %.2f and %.2f" % (MSA_filter_tool, RF_diff_mean, RF_diff_median, RF_diff_std))
    if (sum(df['RF_change'] == "unchanged") + sum(df['RF_change'] == "better")) > 0.8:
      MSA_filter_tools_best.append(MSA_filter_tool)
  
  df_best_filter_tools = df.loc[df["MSA_filter_tool"].isin(MSA_filter_tools_best), ["MSA_filter_tool", "RF_distance_diff", "RF_change"]]
  df_best_filter_tools_automated2 = df_best_filter_tools.loc[df_best_filter_tools["MSA_filter_tool"].str.contains("automated2") , ["MSA_filter_tool", "RF_distance_diff", "RF_change"]]

  #boxplot = sns.boxplot(x=df_best_filter_tools["MSA_filter_tool"], y=df_best_filter_tools["RF_distance_diff"])
  #boxplot.set_xticklabels(boxplot.get_xticklabels(),rotation=90)

  #violinplot= sns.violinplot(x=df_best_filter_tools["MSA_filter_tool"], y=df_best_filter_tools["RF_distance_diff"])
  #violinplot.set_xticklabels(violinplot.get_xticklabels(),rotation=90)

  #joyplot(data=df_best_filter_tools[['RF_distance_diff', 'MSA_filter_tool']], by='MSA_filter_tool',figsize=(12, 8))

  #print(type(df_best_filter_tools_automated2.groupby('MSA_filter_tool')['RF_change'].value_counts(normalize=True) * 100))
  #df_grouped = df_best_filter_tools.groupby('MSA_filter_tool')['RF_change']
  #df_grouped['RF_change_perc'] = df_best_filter_tools.groupby('MSA_filter_tool')['RF_change'].value_counts(normalize=True) * 100
  #df_best_filter_tools['RF_change_perc'] = df_best_filter_tools['RF_change'].value_counts(normalize=True) * 100
  
  #df_best_filter_tools.groupby('MSA_filter_tool')['RF_change'].value_counts(normalize=True).plot.bar()

  #data = df_best_filter_tools.groupby('MSA_filter_tool')['RF_change']

  #bar1 = sns.barplot(x="MSA_filter_tool",  y="RF_change", data=data, color='darkblue')

  #ax = sns.countplot(x="MSA_filter_tool", data=df_best_filter_tools, hue="RF_change")

  fig, ax = plt.subplots(figsize=(8,10))
  #df_best_filter_tools.groupby(['MSA_filter_tool', 'RF_change']).size().unstack().plot(ax=ax, kind='bar', stacked=True, sort_columns=True)
  df_best_filter_tools.groupby(['MSA_filter_tool'])['RF_change'].value_counts(normalize=True).mul(100).unstack()\
    .plot(ax=ax, kind='bar', stacked=True, sort_columns=True, color=['darkgreen', 'darkorange', 'darkred'])

  #his1 = sns.histplot(data=data, x="MSA_filter_tool", hue="RF_change", multiple="stack", stat="percent")
  #ax.set_xticklabels(ax.get_xticklabels(),rotation=90)


  #sns.histplot(df_best_filter_tools, x='MSA_filter_tool', hue='RF_change', weights='RF_change_perc',
  #           multiple='stack', palette='tab20c', shrink=0.8)
  plt.tight_layout()
  plt.show()





if __name__ == "__main__":
  sys.exit(main())
