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

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import OneHotEncoder
from sklearn.tree import DecisionTreeClassifier, export_graphviz


df = pd.DataFrame()


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("--tree", dest = "tree", default = False, action = "store_true", help = "Build decision tree")
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

  table_filename = "table_all_tools_%s_fixed.csv" % args.residueType
  global df
  if os.path.exists(table_filename) and not args.recalculate:
    df = pd.read_csv(table_filename, index_col = 0)
  else:
    data_folder = "alignment_statistics_%s" % args.residueType
    num_sequences = pd.read_table("%s/number_sequences.txt" % data_folder, header = None, names=["num_sequences"], dtype="int")
    num_blocks = pd.read_table("%s/number_blocks.txt" % data_folder, header = None, names=["num_blocks"], dtype="int")
    left_block_column = pd.read_table("%s/left_block_column.txt" % data_folder, header = None, names=["left_block_column"], dtype="int")
    right_block_column = pd.read_table("%s/right_block_column.txt" % data_folder, header = None, names=["right_block_column"], dtype="int")
    num_columns = pd.read_table("%s/number_columns.txt" % data_folder, header = None, names=["num_columns"], dtype="int")
    msa_columns = pd.read_table("%s/number_columns_MSA.txt" % data_folder, header = None, names=["msa_columns"], dtype="int")
    min_columns = pd.read_table("%s/min_columns.txt" % data_folder, header = None, names=["min_columns"])
    max_columns = pd.read_table("%s/max_columns.txt" % data_folder, header = None, names=["max_columns"])
    avg_gaps = pd.read_table("%s/avg_gaps.txt" % data_folder, header = None, names=["avg_gaps"], dtype="float")
    avg_seq_identity = pd.read_table("%s/avg_seq_identity.txt" % data_folder, header = None, names=["avg_seq_identity"])
    RF_distance = pd.read_table("%s/RF_distance.txt" % data_folder, header = None, names=["RF_distance"])
    residue_type = pd.read_table("%s/residue_type.txt" % data_folder, header = None, names=["residue_type"], dtype="string")
    taxon = pd.read_table("%s/taxon.txt" % data_folder, header = None, names=["taxon"], dtype="string")
    problem_num = pd.read_table("%s/problem_num.txt" % data_folder, header = None, names=["problem_num"], dtype="int")
    error = pd.read_table("%s/error_problem.txt" % data_folder, header = None, names=["error"], dtype="string")
    msa_tools = pd.read_table("%s/MSA_tools.txt" % data_folder, header = None, names=["msa_tools"], dtype="string")
    msa_filter_tools = pd.read_table("%s/MSA_filter_tools.txt" % data_folder, header = None, names=["msa_filter_tools"], dtype="string")
    avg_gaps_diff_weighted = pd.read_table("%s/avg_gaps_diff_weighted.txt" % data_folder, header = None, names=["avg_gaps_diff_weighted"], dtype="float")
    avg_seq_identity_diff_weighted = pd.read_table("%s/avg_seq_identity_diff_weighted.txt" % data_folder, header = None, names=["avg_seq_identity_diff_weighted"], dtype="float")
    blocks_diff = pd.read_table("%s/blocks_diff.txt" % data_folder, header = None, names=["blocks_diff"], dtype="int")
    percent_conserved_columns = pd.read_table("%s/percent_conserved_columns.txt" % data_folder, header = None, names=["percent_conserved_columns"], dtype="float")
    removed_columns = pd.read_table("%s/removed_columns.txt" % data_folder, header = None, names=["removed_columns"], dtype="int")
    RF_distance_diff = pd.read_table("%s/RF_distance_diff.txt" % data_folder, header = None, names=["RF_distance_diff"], dtype="int")
    gappy_columns_50 = pd.read_table("%s/gappy_columns_50.txt" % data_folder, header = None, names=["gappy_columns_50"], dtype="float")
    gappy_columns_80 = pd.read_table("%s/gappy_columns_80.txt" % data_folder, header = None, names=["gappy_columns_80"], dtype="float")

    df = pd.concat([num_sequences, num_blocks, left_block_column, right_block_column, num_columns, msa_columns, min_columns, max_columns, avg_seq_identity,
     avg_gaps, RF_distance, residue_type, taxon, problem_num, error, msa_tools, msa_filter_tools, gappy_columns_50, gappy_columns_80,
     removed_columns, percent_conserved_columns, blocks_diff, avg_gaps_diff_weighted, avg_seq_identity_diff_weighted, RF_distance_diff], axis = 1)
    df["has_block"] = df["num_blocks"] > 0
    df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
    df["perc_main_block_size"] = df["main_block_size"] / df["num_columns"]
  
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

  if not os.path.exists(table_filename) or args.recalculate:
    with open(table_filename, 'w') as file:
      file.write(df[df["residue_type"] == args.residueType].to_csv())


  if args.tree:
    run_decision_tree_classifier(df, args.maxDepth, diff = args.compare, ratio =  args.ratio, residue_type = args.residueType, taxon = args.taxon,
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

  df_model = df_model.loc[:, (features + [class_feature])]

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
  class_names = ['worse', 'unchanged', 'better'] if diff else ['different', 'similar']

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
  graph.write_png(tree_filename)


if __name__ == "__main__":
  sys.exit(main())
