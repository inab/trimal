#!/usr/bin/python

#
#  'parse_RF_distances.py'
#
#   Script implemented to parse RF distances data
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
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pySankey.sankey import sankey


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("-t", "--type", dest = "residueType", required = False, type = str, choices = ["AA", "DNA"], help = "Residue type")
  parser.add_argument("--taxon", dest = "taxon", required = False, type = str, choices = ["Bacteria", "Eukaryotes", "Fungi"], help = "Taxon")
  parser.add_argument("--msa_tool", dest = "msaTool", required = False, type = str, choices = ["ClustalW", "ClustalW2", "Mafft", "Prank", "T-Coffee"], help = "MSA tool")
  parser.add_argument("--msa_filter_tool", dest = "msaFilterTool", required = False, type = str, choices = ["trimAl",], help = "MSA filter tool")
  parser.add_argument("--plot", dest = "plotType", required = False, type = str, choices = ["bar_chart", "sankey"], help = "Plot type")


  args = parser.parse_args()

  df = pd.read_csv("RF_distances_AA_automated3_complete.csv")

  problems_to_ignore = pd.unique(df.loc[df["RF_distance"] == -1, "problem_num"].values)
  df = df.loc[~df["problem_num"].isin(problems_to_ignore), :]
  
  df = df.loc[df["MSA_filter_tool"] != 'None', :]
  #df['RF_distance_diff'] = df['RF_distance_diff'].astype(int)
  df.loc[df["RF_distance_diff"] < 0 , 'RF_change'] = "worse"
  df.loc[df["RF_distance_diff"] > 0 , 'RF_change'] = "better"
  df.loc[df["RF_distance_diff"] == 0 , 'RF_change'] = "unchanged"
  df.loc[(df["RF_distance_diff"] == 0) & (df["RF_distance"] == 0), 'RF_change'] = " identical_ref"

  MSA_filter_tools = pd.unique(df['MSA_filter_tool'])

  MSA_filter_tools_best = []
  max_score = 0
  max_score_better = 0
  MSA_filter_tools_scores = dict()


  #df_set = df.loc[~((df["RF_distance_diff"] == 0) & (df["RF_distance"] == 0)), :]
  df_set = df.loc[:, :]
 
  df_set = df_set.set_index(["taxon", "problem_num", "MSA_tool"])

  df_gappyout = df_set.loc[df_set["MSA_filter_tool"] == "trimAl_gappyout", :]
  df_gappyout = df_gappyout.groupby(["taxon", "problem_num", "MSA_tool"])["RF_change"].transform(lambda x: x)

  df_set = df_set.merge(df_gappyout, how='left', on=["taxon", "problem_num", "MSA_tool"], suffixes=('', '_gappyout'))
  
  
  df_set = df_set.reset_index()
  
  df_stats = pd.read_csv("rest_of_tools_alignment_statistics_AA/table_all_tools_AA_fixed.csv")
  df_stats_unfiltered = df_stats.loc[(df_stats['msa_tools'] != 'None') & (df_stats['msa_filter_tools'] == 'None') & ((df_stats['RF_distance'] % 2) == 0), :]
  df_stats = df_stats.loc[(df_stats['msa_tools'] != 'None') & (df_stats['msa_filter_tools'] != 'None') & ((df_stats['RF_distance_diff'] % 2) == 0), :]

  
  df_stats = df_stats.loc[(df_stats["RF_distance"] >= 0) & (df_stats["msa_filter_tools"] == "trimAl_gappyout"), :]
  df_stats = df_stats.set_index(["taxon", "problem_num", "msa_tools"])
  df_set = df_set.merge(df_stats, how='left', left_on=["taxon", "problem_num", "MSA_tool"], right_on=["taxon", "problem_num", "msa_tools"], suffixes=('', '_stats'))


  columns_to_select = ["taxon", "problem_num", "MSA_tool", "number_columns", "MSA_filter_tool", "RF_distance", "RF_distance_diff", "RF_change", "RF_change_gappyout",\
    "num_blocks", "num_columns", "avg_seq_identity", "avg_gaps", "RF_distance_stats", "msa_filter_tools", "removed_columns",\
    "percent_conserved_columns", "avg_gaps_diff_weighted", "avg_seq_identity_diff_weighted", "RF_distance_diff_stats", "main_block_size", "perc_main_block_size", "msa_columns"]

  df_set = df_set.loc[:, columns_to_select]
  df_set = df_set.loc[df_set["problem_num"] != 10, :]
  df_set.loc[:, 'percent_conserved_columns_filtered'] = df_set.loc[:, 'number_columns'] / df_set.loc[:, 'msa_columns'] # here new atuoamted merthods are included
  df_set.loc[:, 'percent_conserved_columns_diff'] = df_set.loc[:, 'percent_conserved_columns'] - df_set.loc[:, 'percent_conserved_columns_filtered']
  #df_set = df_set.loc[(df_set["msa_columns"] > 250), :]
  

  if args.plotType == "bar_chart":
    for MSA_filter_tool in MSA_filter_tools:
      break
      
      #RF_diff_mean = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].mean()
      #RF_diff_median = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].median()
      #RF_diff_std = df.loc[(df['MSA_filter_tool'] == MSA_filter_tool), 'RF_distance_diff'].std()
      #print("Mean, median and std RF_diff for %s are %.2f, %.2f and %.2f" % (MSA_filter_tool, RF_diff_mean, RF_diff_median, RF_diff_std))
      current_score = sum((df_set['MSA_filter_tool'] == MSA_filter_tool) & (df_set['RF_change'] == "unchanged")) /\
        sum(df_set['MSA_filter_tool'] == MSA_filter_tool)

      current_score_better_unchaged_right = (sum((df_set['MSA_filter_tool'] == MSA_filter_tool) & (df_set['RF_change'] == "aidentical_ref")) +\
        sum((df_set['MSA_filter_tool'] == MSA_filter_tool) & (df_set['RF_change'] == "better"))) /\
        sum(df_set['MSA_filter_tool'] == MSA_filter_tool)

      current_score_better = (sum((df_set['MSA_filter_tool'] == MSA_filter_tool) & (df_set['RF_change'] == "better"))) /\
        sum(df_set['MSA_filter_tool'] == MSA_filter_tool)
      print("Perc of better trees of %s = %.2f" % (MSA_filter_tool, current_score_better))
      print("Perc of better and unchanged(right) trees of %s = %.2f" % (MSA_filter_tool, current_score_better_unchaged_right))
      print("Perc of better and unchanged trees of %s = %.2f" % (MSA_filter_tool, current_score+current_score_better_unchaged_right))
      pass
      MSA_filter_tools_scores.setdefault(MSA_filter_tool, current_score)
      if current_score > max_score:
        max_score = current_score
      if current_score_better > max_score_better:
        max_score_better = current_score_better
        #MSA_filter_tools_best.append(MSA_filter_tool)
  
  #df_best_filter_tools = df.loc[df["MSA_filter_tool"].isin(MSA_filter_tools_best), ["MSA_filter_tool", "RF_distance_diff", "RF_change"]]
  #df_best_filter_tools_automated2 = df.loc[df["MSA_filter_tool"].str.contains("automated2") , ["MSA_filter_tool", "RF_distance_diff", "RF_change"]]
  
  '''
  MSA_filter_tools_scores_sorted = dict(sorted(MSA_filter_tools_scores.items(), key=lambda kv: kv[1]))
  
  with open('MSA_filter_tools_scores_sorted_fun.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=MSA_filter_tools_scores_sorted.keys())
    writer.writeheader()
    writer.writerow(MSA_filter_tools_scores_sorted)
  '''


  '''
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
    'trimAl_automated2_similarity_45_250', 'trimAl_automated3']
  '''

  MSA_filter_tools_order = ['Aliscore', 'BMGE62', 'Gblocks', 'Gblocks_relaxed', 'Guidance', 'Noisy08', 'Zorro4', 'trimAl_automated1',\
    'trimAl_gappyout', 'trimAl_strict', 'trimAl_automated2_gaps_30_100',\
    'trimAl_automated2_gaps_30_150', 'trimAl_automated2_gaps_30_200', 'trimAl_automated2_gaps_30_250',\
    'trimAl_automated2_gaps_35_100', 'trimAl_automated2_gaps_35_150', 'trimAl_automated2_gaps_35_200',\
    'trimAl_automated2_gaps_35_250', 'trimAl_automated2_gaps_40_100', 'trimAl_automated2_gaps_40_150',\
    'trimAl_automated2_gaps_40_200', 'trimAl_automated2_gaps_40_250', 'trimAl_automated2_gaps_45_100',\
    'trimAl_automated2_gaps_45_150', 'trimAl_automated2_gaps_45_200', 'trimAl_automated2_gaps_45_250',\
    'trimAl_automated3']


  df_grouped = df_set.groupby(['MSA_filter_tool'])['RF_change'].value_counts(normalize=True).mul(100).unstack()
  df_grouped_sorted = df_grouped.loc[MSA_filter_tools_order]
  
  if args.plotType == "bar_chart":
    fig, ax = plt.subplots(figsize=(17,10))
    #df_best_filter_tools.groupby(['MSA_filter_tool', 'RF_change']).size().unstack().plot(ax=ax, kind='bar', stacked=True, sort_columns=True)
    ticks_x = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    df_grouped_sorted.plot(ax=ax, kind='barh', stacked=True, sort_columns=False, color=['darkblue', 'darkgreen', 'darkorange', 'darkred'], xticks=ticks_x)

    min_score_worse = sum((df_set['MSA_filter_tool'] == 'trimAl_gappyout') & (df_set['RF_change'] == "unchanged")) /\
            sum(df_set['MSA_filter_tool'] == 'trimAl_gappyout')

    max_score_better = (sum((df_set['MSA_filter_tool'] == 'trimAl_gappyout') & (df_set['RF_change'] == " identical_ref")) +\
      sum((df_set['MSA_filter_tool'] == 'trimAl_gappyout') & (df_set['RF_change'] == "better"))) /\
      sum(df_set['MSA_filter_tool'] == 'trimAl_gappyout')

    threshold = (min_score_worse + max_score_better) * 100
    ax.plot([threshold, threshold], [0., len(MSA_filter_tools_order)], "k--")
    threshold = max_score_better * 100
    ax.plot([threshold, threshold], [0., len(MSA_filter_tools_order)], "k--")

    min_score_worse = sum((df_set['MSA_filter_tool'] == 'trimAl_automated3') & (df_set['RF_change'] == "unchanged")) /\
            sum(df_set['MSA_filter_tool'] == 'trimAl_automated3')

    max_score_better = (sum((df_set['MSA_filter_tool'] == 'trimAl_automated3') & (df_set['RF_change'] == " identical_ref")) +\
      sum((df_set['MSA_filter_tool'] == 'trimAl_automated3') & (df_set['RF_change'] == "better"))) /\
      sum(df_set['MSA_filter_tool'] == 'trimAl_automated3')

    threshold = (min_score_worse + max_score_better) * 100
    ax.plot([threshold, threshold], [0., len(MSA_filter_tools_order)], "k--", color='w')
    threshold = max_score_better * 100
    ax.plot([threshold, threshold], [0., len(MSA_filter_tools_order)], "k--", color='w')
    #his1 = sns.histplot(data=data, x="MSA_filter_tool", hue="RF_change", multiple="stack", stat="percent")
    #ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    ax.invert_yaxis()
    ax.set_title('AA')
    ax.set_ylabel('')
    plt.legend(loc='center', title='')
    plt.tight_layout()
    plt.show()
  elif args.plotType == "sankey":
    colors = {
        "worse": "darkred",
        "unchanged": "darkorange",
        "better": "darkgreen"
    }

    df_gappyout = df_set.loc[df_set["MSA_filter_tool"] == "trimAl_gappyout", :]
    df_automated3 = df_set.loc[df_set["MSA_filter_tool"] == "trimAl_automated3", :]
    sankey(df_gappyout["RF_change"], df_automated3["RF_change"], colorDict=colors, aspect=20, fontsize=12,\
       leftLabels=["better", "unchanged", "worse"], rightLabels=["better", "unchanged", "worse"], figure_name="test")

    plt.title("trimAl_gappyout                                              df_automated3")
    plt.tight_layout()
    #plt.savefig("{}.png".format("sankey_diagrams/AA/All/from_1000_columns/trimAl_gappyout-automated2"), bbox_inches='tight', dpi=150)
    plt.show()

    return
    
    df_gappyout = df_set.loc[df_set["MSA_filter_tool"] == "trimAl_gappyout", :]
    df_automated2 = df_set.loc[df_set["MSA_filter_tool"].str.contains("automated2"), :]
    sankey(df_gappyout["RF_change"], df_automated2["RF_change"], colorDict=colors, aspect=20, fontsize=12,\
       leftLabels=["better", "unchanged", "worse"], rightLabels=["better", "unchanged", "worse"], figure_name="test")

    plt.title("trimAl_gappyout                                              automated2(all)")
    plt.tight_layout()
    #plt.savefig("{}.png".format("sankey_diagrams/AA/All/from_1000_columns/trimAl_gappyout-automated2"), bbox_inches='tight', dpi=150)
    plt.show()

    
    for method in filter(lambda method_type: "automated2" in method_type, MSA_filter_tools_order):
      df_automated2_method = df_set.loc[df_set["MSA_filter_tool"] == method, :]
      left = df_gappyout["RF_change"].reset_index(drop=True)
      right = df_automated2_method["RF_change"].reset_index(drop=True)
      df_common = pd.DataFrame({'left': left, 'right': right},\
       index=range(len(df_gappyout["RF_change"])))
      df_common.dropna(inplace=True)
      print(df_common.info())
      
      sankey(df_common.left, df_common.right, colorDict=colors, aspect=20, fontsize=12,\
        leftLabels=["better", "unchanged", "worse"], rightLabels=["better", "unchanged", "worse"])
      
      plt.title("trimAl_gappyout                                              " + method)
      plt.tight_layout()
      #plt.savefig("{}.png".format("sankey_diagrams/AA/All/from_1000_columns/trimAl_gappyout_" + method), bbox_inches='tight', dpi=150)
      plt.show()
  

  
  #df_automated2 = df.loc[df["MSA_filter_tool"] == "trimAl_automated2_gaps_45_250" , :]
  #df_gappyout = df.loc[df["MSA_filter_tool"] == "trimAl_gappyout", :]

  #pandarallel.initialize(progress_bar=True)

  #df_gappyout_values = df.parallel_apply(lambda x: df.loc[(df['taxon'] == x['taxon']) & (df['problem_num'] == x['problem_num']) & (df['MSA_tool'] == x['MSA_tool']) & (df['MSA_filter_tool'] == 'trimAl_gappyout'),\
  #   'RF_change'].values[0], axis=1)

  #df_gappyout_values = df.iloc[:20].apply(lambda x: df.loc[(df['taxon'] == x['taxon']) & (df['problem_num'] == x['problem_num']) & (df['MSA_tool'] == x['MSA_tool']) & (df['MSA_filter_tool'] == 'trimAl_gappyout'),\
  #  'RF_change'].values[0], axis=1)


  else:
    '''
    
    sns.histplot(data=df_set, x="RF_change_gappyout", hue="RF_change", multiple="dodge", shrink=.8, palette=colors, log_scale=(False, True))
    plt.show()

    sns.barplot(data=df_set, x="RF_change_gappyout", y="number_columns", hue="RF_change", palette=colors)
    plt.show()
    '''

    colors = {
        "worse": "darkred",
        "unchanged": "darkorange",
        "better": "darkgreen"
    }


    #print(df_set.head())
    #print(df_set.describe(include='all'))




    df_set = df_set.loc[df_set["MSA_filter_tool"].str.contains("automated3"), :]
    print(df_stats_unfiltered.info())
    #print(df_set.describe())
    print(df_stats_unfiltered["RF_distance"].value_counts(normalize=True))

    sns.histplot(data=df_stats_unfiltered, x="RF_distance", hue="msa_tools", multiple="dodge", stat="percent")
    plt.show()

    print(df_set["RF_distance"].value_counts(normalize=True))
    sns.histplot(data=df_set, x="RF_distance", hue="MSA_filter_tool", multiple="dodge", stat="percent")
    plt.show()

    return

    sns.histplot(data=df_set, x="RF_change_gappyout", hue="RF_change", multiple="dodge", shrink=.8, palette=colors, stat="percent")
    #plt.savefig("{}.png".format("trimAl_automated2_gaps_45_250_plots/AA/All/from_1000_columns/histplot_RF_change"), bbox_inches='tight')
    plt.show()

    sns.violinplot(data=df_set, x="RF_change_gappyout", y="msa_columns", hue="RF_change", palette=colors)
    #plt.savefig("{}.png".format("trimAl_automated2_gaps_45_250_plots/AA/All/from_1000_columns/violinplot_RF_change_msa_columns"), bbox_inches='tight')
    plt.show()

    fig, ax = plt.subplots(ncols=2)
    sns.violinplot(data=df_set, x="RF_change_gappyout", y="percent_conserved_columns", hue="RF_change", palette=colors, ax=ax[0])
    sns.violinplot(data=df_set, x="RF_change_gappyout", y="percent_conserved_columns_filtered", hue="RF_change", palette=colors, ax=ax[1])
    #plt.savefig("{}.png".format("trimAl_automated2_gaps_45_250_plots/AA/All/from_1000_columns/violinplot_RF_change_percent_conserved_columns"), bbox_inches='tight')
    plt.show()

    sns.violinplot(data=df_set, x="RF_change_gappyout", y="percent_conserved_columns_diff", hue="RF_change", palette=colors)
    #plt.savefig("{}.png".format("trimAl_automated2_gaps_45_250_plots/AA/All/from_1000_columns/violinplot_RF_change_percent_conserved_columns_diff"), bbox_inches='tight')
    plt.show()


    #with open("RF_distances_merged.csv", 'w') as file:
    #  file.write(df_set.to_csv())

    '''
    plt.bar(df_automated2['RF_change'], df_automated2['number_columns'])
    plt.show()
    plt.bar(df_gappyout['RF_change'], df_gappyout['number_columns'])
    plt.show()

    sns.histplot(data=df_automated2, x="MSA_tool", hue="RF_change", multiple="stack")
    plt.show()
    sns.histplot(data=df_gappyout, x="MSA_tool", hue="RF_change", multiple="stack")
    plt.show()

    df_common = df.loc[(df["MSA_filter_tool"] == "trimAl_gappyout") | (df["MSA_filter_tool"] == "trimAl_automated2_gaps_45_250"), :]
    sns.catplot(x="MSA_filter_tool", y="number_columns", hue="RF_change", data=df_common, kind="bar")
    plt.show()

      df_automated2_worse = df.loc[(df["MSA_filter_tool"] == "trimAl_automated2_gaps_45_250") & (df["RF_change"] == "worse"), :]
    df_gappyout_better = df.loc[(df["MSA_filter_tool"] == "trimAl_gappyout") & (df["RF_change"] == "better"), :]
    df_better_to_worse = pd.merge(df_gappyout_better, df_automated2_worse, how='inner', on=["taxon", "problem_num", "MSA_tool"])
    print(df_better_to_worse.head())
    #sns.barplot(x=)
    '''

    

    '''
    print(df.head(), flush=True)
    df['RF_change_gappyout'] = "None"
    # add taxon when joined data
    for problem, MSA_tool in zip(df['problem_num'], df['MSA_tool']):
      RF_change_gappyout = df.loc[(df['problem_num'] == problem) & (df['MSA_tool'] == MSA_tool) & \
        (df['MSA_filter_tool'] == "trimAl_gappyout"), 'RF_change'].values
      if RF_change_gappyout.size > 0:
        RF_change_gappyout = RF_change_gappyout[0]
        df.loc[(df['problem_num'] == problem) & (df['MSA_tool'] == MSA_tool), 'RF_change_gappyout'] = RF_change_gappyout
      
      print("problem_num=%i; MSA_tool=%s; RF_change_gappyout=%s" % (int(problem), MSA_tool, RF_change_gappyout), flush=True)

    print(df.head(), flush=True)
    with open("RF_distances_arrays_clean_Bacteria_changes.csv", 'w') as file:
      file.write(df.to_csv())
    '''


if __name__ == "__main__":
  sys.exit(main())
