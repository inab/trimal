#!/usr/bin/python

#
#  'measure_parameters_relevance.py'
#
#   Script implemented to work with trimAl to analyze parameters statistics
#   and decide which are more relevant to the quality of the alignment.
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
import pandas as pd
import pydotplus
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, mean_absolute_error, mean_squared_error, r2_score, classification_report
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.linear_model import LogisticRegression
from sklearn.utils import resample


df = pd.DataFrame()


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("--tree", dest = "tree", default = False, action = "store_true", help = "Build decision tree")
  parser.add_argument("--regression", dest = "regression", default = False, action = "store_true", help = "Build linear regression")
  parser.add_argument("--compare", dest = "compare", default = False, action = "store_true", help = "Compare filtered results with original MSA")
  parser.add_argument("--recalculate", dest = "recalculate", default = False, action = "store_true", help = "Recalculate values")
  parser.add_argument("--max_depth", dest = "maxDepth", default = 2, type = int, help = "Tree max depth")
  parser.add_argument("--ratio", dest = "ratio", default = False,  action = "store_true", help = "Add ratio parameters to the model")
  parser.add_argument("-t", "--type", dest = "residueType", required = True, type = str, choices = ["AA", "DNA"], help = "Residue type")
  parser.add_argument("--taxon", dest = "taxon", required = False, type = str, choices = ["Bacteria", "Eukaryotes", "Fungi"], help = "Taxon")
  parser.add_argument("--msa_tool", dest = "msaTool", required = False, type = str, choices = ["ClustalW", "ClustalW2", "Mafft", "Prank", "T-Coffee"], help = "MSA tool")
  parser.add_argument("-c", "--criterion", dest = "criterion", required = False, default = "gini", type = str, choices = ["gini", "entropy"], help = "The function to measure the quality of a split")
  parser.add_argument("--min_columns", dest = "minColumns", required = False, type = int,  help = "Minimin number of columns of the MSA alignment (unfiltered)")
  parser.add_argument("--min_seqs", dest = "minSeqs", required = False, type = int,  help = "Minimin number of sequences of the alignment")


  args = parser.parse_args()

  table_filename = "rest_of_tools_alignment_statistics_AA/table_all_tools_AA_fixed.csv"
  global df
  df = pd.read_csv(table_filename)

  
  if args.recalculate and args.ratio:
    df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
    df.loc[df["main_block_size"] < 0, "main_block_size"] = 0
    df["columns/sequence"] = df["num_columns"] / df["num_sequences"]
    df["blocks/columns"] = df["num_blocks"] / df["num_columns"]
    df["perc_main_block_size"] = df["main_block_size"] / df["num_columns"]

  '''
  with open('table.html', 'w') as file:
    file.write(df.to_html())

  '''

  df["has_block"] = df["num_blocks"] > 0
  df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
  df["perc_main_block_size"] = df["main_block_size"] / df["num_columns"]
  df = df.loc[(df['msa_tools'] != 'None') & (df['msa_filter_tools'] != 'None') & ((df['RF_distance_diff'] % 2) == 0), :]

  df_better = df[df["RF_distance_diff"] > 0]
  df_worse = df[df["RF_distance_diff"] < 0]
  df_unchanged = df[df["RF_distance_diff"] == 0]
  df_plot = df.loc[(df['msa_tools'] != 'None') & (df['msa_filter_tools'] != 'None') & (df['RF_distance'] != -1), :]

  colors = {
    "worse": "darkred",
    "unchanged": "darkorange",
    "better": "darkgreen"
  }

  df_plot.loc[:, "RF_change"] = "None"
  df_plot.loc[(df_plot['RF_distance_diff'] > 0), 'RF_change'] = "better"
  df_plot.loc[(df_plot['RF_distance_diff'] == 0), 'RF_change'] = "unchanged"
  df_plot.loc[(df_plot['RF_distance_diff'] < 0), 'RF_change'] = "worse"

  print(stats.ttest_ind(a=df_plot.loc[(df_plot['RF_distance_diff'] > 0), 'percent_conserved_columns'], b=df_plot.loc[(df_plot['RF_distance_diff'] == 0), 'percent_conserved_columns'],\
     equal_var=True))
  print(stats.ttest_ind(a=df_plot.loc[(df_plot['RF_distance_diff'] > 0), 'percent_conserved_columns'], b=df_plot.loc[(df_plot['RF_distance_diff'] < 0), 'percent_conserved_columns'],\
     equal_var=True))

  print(stats.ttest_ind(a=df_plot.loc[(df_plot['RF_distance_diff'] == 0), 'percent_conserved_columns'], b=df_plot.loc[(df_plot['RF_distance_diff'] < 0), 'percent_conserved_columns'],\
      equal_var=True))
  
  sns.violinplot(data=df_plot, x="RF_change", y="percent_conserved_columns", palette=colors)
  plt.xlabel("RF change")
  plt.ylabel("% of conserved columns")
  plt.show()

  sns.boxplot(data=df_plot, x="RF_change", y="percent_conserved_columns", palette=colors)
  plt.xlabel("RF change")
  plt.ylabel("% of conserved columns")
  plt.show()

  df_plot = df_plot.loc[df_plot["num_columns"] < 1000, :]


  sns.violinplot(data=df_plot, x="RF_change", y="num_columns", palette=colors)
  plt.xlabel("RF change")
  plt.ylabel("Columns")
  plt.show()

  sns.boxplot(data=df_plot, x="RF_change", y="num_columns", palette=colors)
  plt.xlabel("RF change")
  plt.ylabel("Columns")
  plt.show()
  




  if args.tree:
    run_decision_tree_classifier(df, args.maxDepth, diff = args.compare, ratio =  args.ratio, residue_type = args.residueType, taxon = args.taxon,
    tool = args.msaTool, criterion = args.criterion, min_columns = args.minColumns, min_seqs = args.minSeqs)
  if args.regression:
    log_regression(df, args.maxDepth, diff = args.compare, ratio =  args.ratio, residue_type = args.residueType, taxon = args.taxon,
    tool = args.msaTool, criterion = args.criterion, min_columns = args.minColumns, min_seqs = args.minSeqs)


def run_decision_tree_classifier(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
  features = ['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'main_block_size', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted', \
    "avg_seq_identity_diff_weighted"]
  class_feature = 'RF_distance_diff' if diff else 'RF_distance'
  if diff:
    features += [ 'percent_conserved_columns', 'removed_columns']
  if ratio:
    features += ['columns/sequence', 'blocks/columns']
  df_model = df.copy()
  if taxon:
    df_model = df_model[df_model["taxon"] == taxon]
  if tool:
    df_model = df_model.loc[(df_model['msa_tools'] == tool) & (df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
  else:
    df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
  if min_columns:
    df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
  if min_seqs:
    df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
  if diff:
    df_model.loc[(df_model['RF_distance_diff'] > 0), 'RF_distance_diff'] = 2
    df_model.loc[(df_model['RF_distance_diff'] == 0), 'RF_distance_diff'] = 1
    df_model.loc[(df_model['RF_distance_diff'] < 0), 'RF_distance_diff'] = 0
    df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 3
  else:
    df_model['RF_distance'] = df_model['RF_distance'] < 4
  
  enc = OneHotEncoder()
  enc_df = pd.DataFrame(enc.fit_transform(df_model[['msa_filter_tools']]).toarray())
  enc_df.columns = enc.get_feature_names_out()
  df_model = df_model.reset_index(drop=True)
  df_model = df_model.join(enc_df)
  features += enc.get_feature_names_out().tolist()

  if not tool:
    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(df_model[['msa_tools']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()

  print(df_model.info())
  df_model = df_model.loc[:, (features + [class_feature])]
  df_model = df_model.dropna()

  print(df_model[features].head())
  print(df_model[features].describe())
  print(df_model[features].info())

  print(df_model[class_feature].info())

  model = DecisionTreeClassifier(max_depth = max_depth, criterion = criterion)
  X = df_model[features]
  y = df_model[class_feature]
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle = True, train_size = 0.8)
  model.fit(X_train, y_train)
  print(len(X_train))
  print(len(X_test))

  predictions = model.predict(X_test)

  print("Accuracy:", accuracy_score(y_test, predictions))
  class_names = ['worse', 'unchanged-wrong', 'better', 'unchanged-right'] if diff else ['different', 'similar']

  dot_data = export_graphviz(model, filled = True, rounded = True, special_characters = True, proportion = True, precision = 2,
    feature_names = df_model[features].columns, class_names = class_names)
  graph = pydotplus.graph_from_dot_data(dot_data)
  tree_filename_prex = "tree_%s_" % residue_type
  if taxon:
    tree_filename_prex += "%s_" % taxon
  if tool:
    tree_filename_prex += "%s_" % tool
  if diff:
    tree_filename_prex += "diff_"
  if ratio:
    tree_filename_prex += "ratio_"
  if min_columns:
    tree_filename_prex += "%s_min_columns_" % min_columns
  if min_seqs:
    tree_filename_prex += "%s_min_seqs_" % min_seqs
  tree_filename = tree_filename_prex + str(max_depth) + "_" + criterion + ".png"
  print("saved " + tree_filename)
  graph.write_png(tree_filename)



def log_regression(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
  features = ['num_sequences', 'msa_columns', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'perc_main_block_size', 'avg_gaps_diff_weighted', \
    "avg_seq_identity_diff_weighted"]
  class_feature = 'RF_distance_diff' if diff else 'RF_distance'
  if diff:
    features += [ 'percent_conserved_columns', 'removed_columns']
  if ratio:
    features += ['columns/sequence', 'blocks/columns']
  df_model = df.copy()
  if taxon:
    df_model = df_model[df_model["taxon"] == taxon]
  if tool:
    df_model = df_model.loc[(df_model['msa_tools'] == tool) & (df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
  else:
    df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
  if min_columns:
    df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
  if min_seqs:
    df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
  if diff:
    df_model.loc[(df_model['RF_distance_diff'] > 0), 'RF_distance_diff'] = 2
    df_model.loc[(df_model['RF_distance_diff'] == 0), 'RF_distance_diff'] = 1
    df_model.loc[(df_model['RF_distance_diff'] < 0), 'RF_distance_diff'] = 0
  else:
    df_model['RF_distance'] = df_model['RF_distance'] < 4

  '''
  enc = OneHotEncoder()
  enc_df = pd.DataFrame(enc.fit_transform(df_model[['MSA_filter_tool']]).toarray())
  enc_df.columns = enc.get_feature_names_out()
  df_model = df_model.reset_index(drop=True)
  df_model = df_model.join(enc_df)
  features += enc.get_feature_names_out().tolist()

  if not tool:
    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(df_model[['MSA_tool']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()



  '''


  df_model = df_model.dropna()
  std_scaler = StandardScaler()
  df_model = df_model.loc[:, (features + [class_feature])]
  df_model.loc[:, features] = std_scaler.fit_transform(df_model.loc[:, features])
  


  print(df_model[features].head())
  print(df_model[features].describe())
  print(df_model[features].info())

  print(df_model[class_feature].info())

  model = LogisticRegression()
  X = df_model[features]
  y = df_model[class_feature]
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle = True, train_size = 0.8)
  model.fit(X_train, y_train)
  print(len(X_train))
  print(len(X_test))

  print(model.score(X_train, y_train))

  cdf = pd.concat([pd.DataFrame(X.columns),pd.DataFrame(np.transpose(model.coef_))], axis = 1)
  print(cdf)

  predictions = model.predict(X_test)

  mae = mean_absolute_error(y_test, predictions)
  mse = mean_squared_error(y_test, predictions)
  r2 = r2_score(y_test, predictions)

  print(classification_report(y_test,predictions))

  print("The model performance for testing set")
  print("--------------------------------------")
  print('MAE is {}'.format(mae))
  print('MSE is {}'.format(mse))
  print('R2 score is {}'.format(r2))

  corrMatrix = df_model.corr()
  sns.heatmap(corrMatrix, annot=True)
  plt.show()
  #print(corrMatrix)


if __name__ == "__main__":
  sys.exit(main())
