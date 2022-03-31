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
import numpy as np
import pandas as pd
import pydotplus

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from joblib import Parallel, delayed


df = pd.DataFrame()


def main():
  parser = argparse.ArgumentParser()

  parser.add_argument("--compare_prank", dest = "comparePrank", default = False, action = "store_true", help = "Compare filtered results with Prank's")
  parser.add_argument("--recalculate", dest = "recalculate", default = False, action = "store_true", help = "Recalculate values")
  parser.add_argument("--max_depth", dest = "maxDepth", default = 2, type = int, help = "Tree max depth")
  parser.add_argument("--ratio", dest = "ratio", default = False,  action = "store_true", help = "Add ratio parameters to the model")
  parser.add_argument("-t", "--type", dest = "residueType", required = True, type = str, choices = ["AA", "DNA"], help = "Residue type")
  parser.add_argument("--taxon", dest = "taxon", required = False, type = str, choices = ["Bacteria", "Eukaryotes", "Fungi"], help = "Taxon")

  args = parser.parse_args()

  table_filename = "table_%s.csv" % args.residueType
  global df
  if os.path.exists(table_filename) and not args.recalculate:
    df = pd.read_csv(table_filename, index_col = 0)
  else:
    num_sequences = pd.read_table("test_working_files/number_sequences.txt", header = None, names=["num_sequences"], dtype="int")
    num_blocks = pd.read_table("test_working_files/blocks_outputs.txt", header = None, names=["num_blocks"], dtype="int")
    left_block_column = pd.read_table("test_working_files/left_block_column.txt", header = None, names=["left_block_column"], dtype="int")
    right_block_column = pd.read_table("test_working_files/right_block_column.txt", header = None, names=["right_block_column"], dtype="int")
    num_columns = pd.read_table("test_working_files/number_columns.txt", header = None, names=["num_columns"], dtype="int")
    min_columns = pd.read_table("test_working_files/min_columns.txt", header = None, names=["min_columns"])
    max_columns = pd.read_table("test_working_files/max_columns.txt", header = None, names=["max_columns"])
    avg_gaps = pd.read_table("test_working_files/avg_gaps.txt", header = None, names=["avg_gaps"], dtype="float")
    avg_seq_identity = pd.read_table("test_working_files/avg_seq_identity.txt", header = None, names=["avg_seq_identity"])
    RF_distance = pd.read_table("test_working_files/RF_distance.txt", header = None, names=["RF_distance"])
    residue_type = pd.read_table("test_working_files/residue_type.txt", header = None, names=["residue_type"], dtype="string")
    taxon = pd.read_table("test_working_files/taxon.txt", header = None, names=["taxon"], dtype="string")
    problem_num = pd.read_table("test_working_files/problem_num.txt", header = None, names=["problem_num"], dtype="int")
    error = pd.read_table("test_working_files/error_problem.txt", header = None, names=["error"], dtype="string")
    msa_tools = pd.read_table("test_working_files/MSA_tools.txt", header = None, names=["msa_tools"], dtype="string")
    msa_filter_tools = pd.read_table("test_working_files/MSA_filter_tools.txt", header = None, names=["msa_filter_tools"], dtype="string")

    df = pd.concat([num_sequences, num_blocks, left_block_column, right_block_column, num_columns, min_columns, max_columns, avg_seq_identity, avg_gaps, RF_distance, residue_type, taxon, problem_num, error, msa_tools, msa_filter_tools], axis = 1)
    df["has_block"] = df["num_blocks"] > 0
  
  if args.recalculate:
    df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
    df.loc[df["main_block_size"] < 0, "main_block_size"] = 0
    df["columns/sequence"] = df["num_columns"] / df["num_sequences"]
    df["blocks/columns"] = df["num_blocks"] / df["num_columns"]
    df["perc_main_block_size"] = df["main_block_size"] / df["num_columns"]
  if args.comparePrank and args.recalculate:
    df["removed_columns"] = -1
    df["perc_conserved_columns"] = -1
    df["blocks_diff"] = np.NaN
    df["avg_gaps_diff"] = np.NaN
    df["avg_gaps_diff_weighted"] = np.NaN
    df["avg_seq_identity_diff"] = np.NaN
    df["RF_distance_diff"] = np.NaN
    #df["avg_seq_identity_diff_weighted"] = df["avg_seq_identity"] * df["perc_conserved_columns"]
    Parallel(n_jobs = 8, require='sharedmem')(delayed(add_comparisons_with_prank)(problem) for problem in range(1, 1000))

  '''
  with open('table.html', 'w') as file:
    file.write(df.to_html())

  '''

  if not os.path.exists(table_filename) or args.recalculate:
    with open(table_filename, 'w') as file:
      file.write(df[df["residue_type"] == args.residueType].to_csv())


  run_decision_tree_classifier(df, args.maxDepth, diff = args.comparePrank, ratio =  args.ratio, residue_type = args.residueType, taxon = args.taxon)


def add_comparisons_with_prank(problem):
  msa_filter_tools_unique = set(df["msa_filter_tools"])
  df_prank = df[df['msa_tools'] == 'Prank']
  df_problem = df_prank[df_prank['problem_num'] == problem]
  for residue in ['AA', 'DNA']:
    df_residue = df_problem[df_problem['residue_type'] == residue]
    for taxa in ['Bacteria', 'Eukaryotes', 'Fungi']:
      df_taxa = df_residue[df_residue['taxon'] == taxa]
      prank_columns = df_taxa.loc[(df_taxa['msa_filter_tools'] == 'None'), ['num_columns']]
      prank_blocks = df_taxa.loc[(df_taxa['msa_filter_tools'] == 'None'), ['num_blocks']]
      prank_avg_gaps = df_taxa.loc[(df_taxa['msa_filter_tools'] == 'None'), ['avg_gaps']]
      prank_avg_seq_identity = df_taxa.loc[(df_taxa['msa_filter_tools'] == 'None'), ['avg_seq_identity']]
      prank_avg_RF_distance = df_taxa.loc[(df_taxa['msa_filter_tools'] == 'None'), ['RF_distance']]
      if prank_columns.values.size > 0:
        for filter_tool in msa_filter_tools_unique:
          if filter_tool != 'None' and filter_tool != 'Original':
            filter_tool_columns = df_taxa.loc[((df_taxa['msa_filter_tools'] == filter_tool)), ['num_columns']]
            filter_tool_blocks = df_taxa.loc[((df_taxa['msa_filter_tools'] == filter_tool)), ['num_blocks']]
            filter_tool_avg_gaps = df_taxa.loc[((df_taxa['msa_filter_tools'] == filter_tool)), ['avg_gaps']]
            filter_tool_avg_seq_identity = df_taxa.loc[((df_taxa['msa_filter_tools'] == filter_tool)), ['avg_seq_identity']]
            filter_tool_RF_distance = df_taxa.loc[((df_taxa['msa_filter_tools'] == filter_tool)), ['RF_distance']]
            if filter_tool_columns.values.size > 0:
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['removed_columns']] = prank_columns.values[0][0] - filter_tool_columns.values[0][0]
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['perc_conserved_columns']] = filter_tool_columns.values[0][0] / prank_columns.values[0][0]
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['blocks_diff']] = prank_blocks.values[0][0] - filter_tool_blocks.values[0][0]
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['avg_gaps_diff']] = prank_avg_gaps.values[0][0] - filter_tool_avg_gaps.values[0][0]
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['avg_gaps_diff_weighted']] = (prank_avg_gaps.values[0][0] - filter_tool_avg_gaps.values[0][0]) \
                * (filter_tool_columns.values[0][0] / prank_columns.values[0][0])
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['avg_seq_identity_diff']] = prank_avg_seq_identity.values[0][0] - filter_tool_avg_seq_identity.values[0][0]
              df.loc[((df['msa_tools'] == 'Prank') & (df['taxon'] == taxa) & (df['residue_type'] == residue) & (df['msa_filter_tools'] == filter_tool) &
                (df['problem_num'] == problem)), ['RF_distance_diff']] = prank_avg_RF_distance.values[0][0] - filter_tool_RF_distance.values[0][0]
              print('Processed problem ' + str(problem) + '-' + residue + '-' + taxa + '-' + filter_tool)


def run_decision_tree_classifier(df, max_depth, diff, ratio, residue_type, taxon):
  features = ['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'main_block_size', 'msa_filter_tools']
  class_feature = 'RF_distance_diff' if diff else 'RF_distance'
  if diff:
    features += ['blocks_diff', 'avg_gaps_diff_weighted', 'avg_seq_identity_diff_weighted', 'perc_conserved_columns', 'removed_columns']
  if ratio:
    features += ['columns/sequence', 'blocks/columns', 'perc_main_block_size']
  df_model = df.loc[:, (features + [class_feature])]
  if taxon:
    df_model = df_model[df["taxon"] == taxon]
  df_model = df_model.loc[(df['msa_tools'] == 'Prank') & (df['msa_filter_tools'] != 'None') & (df['RF_distance'] != -1), :]
  if diff:
    df_model.loc[(df_model['RF_distance_diff'] > 0), 'RF_distance_diff'] = 2
    df_model.loc[(df_model['RF_distance_diff'] < 0), 'RF_distance_diff'] = 1
  else:
    df_model['RF_distance'] = df_model['RF_distance'] <= 4
  
  enc = OneHotEncoder()
  enc_df = pd.DataFrame(enc.fit_transform(df_model[['msa_filter_tools']]).toarray())
  enc_df.columns = enc.get_feature_names_out()
  df_model = df_model.reset_index(drop=True)
  df_model = df_model.join(enc_df)
  features.remove('msa_filter_tools')
  features += enc.get_feature_names_out().tolist()

  print(df_model[features].head())
  print(df_model[features].describe())
  print(df_model[features].info())

  model = DecisionTreeClassifier(max_depth = max_depth)
  X = df_model[features]
  y = df_model[class_feature]
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle = True, train_size = 0.75)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)

  print("Accuracy:", accuracy_score(y_test, predictions))
  class_names = ['equal', 'worse', 'better'] if diff else ['different', 'similar']

  dot_data = export_graphviz(model, filled = True, rounded = True, special_characters = True, proportion = True, precision = 2,
    feature_names = df_model[features].columns, class_names = class_names)
  graph = pydotplus.graph_from_dot_data(dot_data)
  tree_filename_prex = "tree_%s_" % residue_type
  if taxon:
    tree_filename_prex += "%s_" % taxon
  if diff:
    tree_filename_prex += "diff_"
  if ratio:
    tree_filename_prex += "ratio_"
  tree_filename = tree_filename_prex + str(max_depth) + ".png"
  graph.write_png(tree_filename)


if __name__ == "__main__":
  sys.exit(main())
