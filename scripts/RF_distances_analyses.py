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
import csv


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("-t", "--type", dest = "residueType", required = False, type = str, choices = ["AA", "DNA"], help = "Residue type")
  parser.add_argument("--taxon", dest = "taxon", required = False, type = str, choices = ["Bacteria", "Eukaryotes", "Fungi"], help = "Taxon")
  parser.add_argument("--msa_tool", dest = "msaTool", required = False, type = str, choices = ["ClustalW", "ClustalW2", "Mafft", "Prank", "T-Coffee"], help = "MSA tool")
  parser.add_argument("--msa_filter_tool", dest = "msaFilterTool", required = False, type = str, choices = ["trimAl",], help = "MSA filter tool")


  args = parser.parse_args()

  df = pd.read_csv("RF_distances_arrays_fun_clean.csv", index_col = 0)
  df = df.loc[df["RF_distance"] != -1, :]
  df = df.loc[df["MSA_filter_tool"] != 'None', :]
  #df['RF_distance_diff'] = df['RF_distance_diff'].astype(int)
  df.loc[df["RF_distance_diff"] < 0 , 'RF_change'] = "worse"
  df.loc[df["RF_distance_diff"] > 0 , 'RF_change'] = "better"
  df.loc[df["RF_distance_diff"] == 0 , 'RF_change'] = "unchanged"

  MSA_filter_tools = pd.unique(df['MSA_filter_tool'])
  MSA_filter_tools_best = []
  max_score = 0
  max_score_better = 0
  MSA_filter_tools_scores = dict()

  for MSA_filter_tool in MSA_filter_tools:
    #break
    #RF_diff_mean = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].mean()
    #RF_diff_median = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].median()
    #RF_diff_std = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].std()
    #print("Mean, median and std RF_diff for %s are %.2f, %.2f and %.2f" % (MSA_filter_tool, RF_diff_mean, RF_diff_median, RF_diff_std))
    current_score = (sum((df['MSA_filter_tool'] == MSA_filter_tool) & (df['RF_change'] == "unchanged")) +\
      sum((df['MSA_filter_tool'] == MSA_filter_tool) & (df['RF_change'] == "better"))) /\
      sum(df['MSA_filter_tool'] == MSA_filter_tool)

    current_score_better = (sum((df['MSA_filter_tool'] == MSA_filter_tool) & (df['RF_change'] == "better"))) /\
      sum(df['MSA_filter_tool'] == MSA_filter_tool)
    #print("Perc of better and unchanged trees of %s = %.2f" % (MSA_filter_tool, current_score))
    MSA_filter_tools_scores.setdefault(MSA_filter_tool, current_score)
    if current_score > max_score:
      max_score = current_score
    if current_score_better > max_score_better:
      max_score_better = current_score_better
      #MSA_filter_tools_best.append(MSA_filter_tool)
  
  #df_best_filter_tools = df.loc[df["MSA_filter_tool"].isin(MSA_filter_tools_best), ["MSA_filter_tool", "RF_distance_diff", "RF_change"]]
  #df_best_filter_tools_automated2 = df.loc[df["MSA_filter_tool"].str.contains("automated2") , ["MSA_filter_tool", "RF_distance_diff", "RF_change"]]
  
  MSA_filter_tools_scores_sorted = dict(sorted(MSA_filter_tools_scores.items(), key=lambda kv: kv[1]))
  #print(MSA_filter_tools_scores_sorted)
  
  with open('MSA_filter_tools_scores_sorted_fun.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=MSA_filter_tools_scores_sorted.keys())
    writer.writeheader()
    writer.writerow(MSA_filter_tools_scores_sorted)
  


  fig, ax = plt.subplots(figsize=(17,10))
  MSA_filter_tools_order=['Aliscore', 'BMGE62', 'Gblocks', 'Gblocks_relaxed', 'Guidance', 'Noisy08', 'Zorro4', 'trimAl_automated1',\
    'trimAl_gappyout', 'trimAl_strict', 'trimAl_automated2_combined_30_100', 'trimAl_automated2_combined_30_150',\
    'trimAl_automated2_combined_30_200', 'trimAl_automated2_combined_30_250', 'trimAl_automated2_combined_35_100',\
    'trimAl_automated2_combined_35_150', 'trimAl_automated2_combined_35_200', 'trimAl_automated2_combined_35_250',\
    'trimAl_automated2_combined_40_100', 'trimAl_automated2_combined_40_150', 'trimAl_automated2_combined_40_200',\
    'trimAl_automated2_combined_40_250', 'trimAl_automated2_combined_45_100', 'trimAl_automated2_combined_45_150',\
    'trimAl_automated2_combined_45_200', 'trimAl_automated2_combined_45_250', 'trimAl_automated2_gaps_30_100',\
    'trimAl_automated2_gaps_30_150', 'trimAl_automated2_gaps_30_200', 'trimAl_automated2_gaps_30_250',\
    'trimAl_automated2_gaps_35_100', 'trimAl_automated2_gaps_35_150', 'trimAl_automated2_gaps_35_200',\
    'trimAl_automated2_gaps_35_250', 'trimAl_automated2_gaps_40_100', 'trimAl_automated2_gaps_40_150',\
    'trimAl_automated2_gaps_40_200', 'trimAl_automated2_gaps_40_250', 'trimAl_automated2_gaps_45_100',\
    'trimAl_automated2_gaps_45_150', 'trimAl_automated2_gaps_45_200', 'trimAl_automated2_gaps_45_250',\
    'trimAl_automated2_similarity_30_100', 'trimAl_automated2_similarity_30_150', 'trimAl_automated2_similarity_30_200',\
    'trimAl_automated2_similarity_30_250', 'trimAl_automated2_similarity_35_100', 'trimAl_automated2_similarity_35_150',\
    'trimAl_automated2_similarity_35_200', 'trimAl_automated2_similarity_35_250', 'trimAl_automated2_similarity_40_100',\
    'trimAl_automated2_similarity_40_150', 'trimAl_automated2_similarity_40_200', 'trimAl_automated2_similarity_40_250',\
    'trimAl_automated2_similarity_45_100', 'trimAl_automated2_similarity_45_150', 'trimAl_automated2_similarity_45_200',\
    'trimAl_automated2_similarity_45_250']

  df_grouped = df.groupby(['MSA_filter_tool'])['RF_change'].value_counts(normalize=True).mul(100).unstack()
  df_grouped_sorted = df_grouped.loc[MSA_filter_tools_order]
  
  #df_best_filter_tools.groupby(['MSA_filter_tool', 'RF_change']).size().unstack().plot(ax=ax, kind='bar', stacked=True, sort_columns=True)
  df_grouped_sorted.plot(ax=ax, kind='barh', stacked=True, sort_columns=False, color=['darkgreen', 'darkorange', 'darkred'])

  threshold = max_score * 100
  ax.plot([threshold, threshold], [0., len(MSA_filter_tools_order)], "k--")
  threshold = max_score_better * 100
  ax.plot([threshold, threshold], [0., len(MSA_filter_tools_order)], "k--")
  #his1 = sns.histplot(data=data, x="MSA_filter_tool", hue="RF_change", multiple="stack", stat="percent")
  #ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
  ax.invert_yaxis()
  ax.set_title('Bacteria')
  ax.set_ylabel('')
  plt.legend(loc='center', title='')



  plt.tight_layout()
  plt.show()





if __name__ == "__main__":
  sys.exit(main())
