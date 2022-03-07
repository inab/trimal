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

import sys
import os
import pandas as pd
import pydotplus

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, explained_variance_score, mean_squared_error, r2_score
from sklearn.preprocessing import MinMaxScaler, OneHotEncoder
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier, export_graphviz


def main():
  if os.path.exists("table.csv"):
    df = pd.read_csv("table.csv", index_col=0)
  else:
    num_sequences = pd.read_table("test_working_files/number_sequences.txt", header=None, names=["num_sequences"], dtype="int")
    num_blocks = pd.read_table("test_working_files/blocks_outputs.txt", header=None, names=["num_blocks"], dtype="int")
    left_block_column = pd.read_table("test_working_files/left_block_column.txt", header=None, names=["left_block_column"], dtype="int")
    right_block_column = pd.read_table("test_working_files/right_block_column.txt", header=None, names=["right_block_column"], dtype="int")
    num_columns = pd.read_table("test_working_files/number_columns.txt", header=None, names=["num_columns"], dtype="int")
    min_columns = pd.read_table("test_working_files/min_columns.txt", header=None, names=["min_columns"])
    max_columns = pd.read_table("test_working_files/max_columns.txt", header=None, names=["max_columns"])
    avg_gaps = pd.read_table("test_working_files/avg_gaps.txt", header=None, names=["avg_gaps"], dtype="float")
    avg_seq_identity = pd.read_table("test_working_files/avg_seq_identity.txt", header=None, names=["avg_seq_identity"])
    RF_distance = pd.read_table("test_working_files/RF_distance.txt", header=None, names=["RF_distance"])
    residue_type = pd.read_table("test_working_files/residue_type.txt", header=None, names=["residue_type"], dtype="string")
    taxon = pd.read_table("test_working_files/taxon.txt", header=None, names=["taxon"], dtype="string")
    problem_num = pd.read_table("test_working_files/problem_num.txt", header=None, names=["problem_num"], dtype="string")
    error = pd.read_table("test_working_files/error_problem.txt", header=None, names=["error"], dtype="string")
    msa_tools = pd.read_table("test_working_files/MSA_tools.txt", header=None, names=["msa_tools"], dtype="string")
    msa_filter_tools = pd.read_table("test_working_files/MSA_filter_tools.txt", header=None, names=["msa_filter_tools"], dtype="string")

    df = pd.concat([num_sequences, num_blocks, left_block_column, right_block_column, num_columns, min_columns, max_columns, avg_seq_identity, avg_gaps, RF_distance, residue_type, taxon, problem_num, error, msa_tools, msa_filter_tools], axis=1)
    df["has_block"] = df["num_blocks"] > 0
  
  df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
  df[df["main_block_size"] < 0] = 0

  print(df.head())
  print(df.tail())
  print(df.info())

  '''
  with open('table.csv', 'w') as file:
    file.write(df.to_csv())
  print(df.corr().to_string())

  with open('table.html', 'w') as file:
    file.write(df.to_html())

  '''
  if not os.path.exists("table.csv"):
    with open('table.csv', 'w') as file:
      file.write(df.to_csv())
  

  run_decision_tree_classifier(df)


def run_regression(df):
  scaler = MinMaxScaler()
  df_model = df[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'RF_distance', 'has_block', 'main_block_size']]
  df_model = df_model[(df['msa_tools'] != 'Original') & (df['RF_distance'] != -1)]
  df_model['is_AA'] = df['residue_type'] == "AA"
  df_model = pd.DataFrame(scaler.fit_transform(df_model),
            columns=df_model.columns, index=df_model.index)
  
  print(df_model.corr().to_string())

  model = LinearRegression()
  X = df_model[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'has_block', 'main_block_size', 'is_AA']]
  y = df_model['RF_distance']
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, train_size=0.7)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)
  r2 = r2_score(y_test, predictions)
  rmse = mean_squared_error(y_test, predictions, squared=False)

  print("Explained variance:", explained_variance_score(y_test, predictions))
  print("r2 score:", r2_score(y_test, predictions))
  print('The r2 is: ', r2)
  print('The rmse is: ', rmse)
  print(model.coef_)
  print(model.intercept_)


def run_decision_tree_regressor(df):
  scaler = MinMaxScaler()
  df_model = df[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'RF_distance', 'has_block', 'main_block_size']]
  df_model = df_model[(df['msa_tools'] != 'Original') & (df['RF_distance'] != -1)]
  df_model['is_AA'] = df['residue_type'] == "AA"
  
  #print(df_model.corr().to_string())
  print(df_model.describe())

  df_model = pd.DataFrame(scaler.fit_transform(df_model),
            columns=df_model.columns, index=df_model.index) 

  print(df_model.describe())

  df_model_features = df_model[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'has_block', 'main_block_size', 'is_AA']]

  model = DecisionTreeRegressor()
  X = df_model_features
  y = df_model['RF_distance']
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, train_size=0.7)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)

  print("Explained variance:", explained_variance_score(y_test, predictions))
  print("r2 score:", r2_score(y_test, predictions))


def run_decision_tree_classifier(df):
  #scaler = MinMaxScaler()
  df_model = df[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'RF_distance', 'has_block', 'main_block_size']]
  df_model = df_model[(df['msa_tools'] != 'Original') & (df['RF_distance'] != -1)]
  df_model['RF_distance'] = df_model['RF_distance'] > 4
  df_model['is_AA'] = df['residue_type'] == "AA"

  df_model_features = df_model[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'has_block', 'main_block_size', 'is_AA']]
  print(df_model_features.describe())
  print(df_model_features.info())
  #print(df_model.corr().to_string())

  model = DecisionTreeClassifier(max_depth=3)
  X = df_model_features
  y = df_model['RF_distance']
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, train_size=0.7)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)

  print("Accuracy:", accuracy_score(y_test, predictions))

  dot_data = export_graphviz(model, filled=True, rounded=True, special_characters=True, feature_names = df_model_features.columns, class_names=['Bad', 'Good'])
  graph = pydotplus.graph_from_dot_data(dot_data)
  graph.write_png('tree.png')


def run_random_forest_classifier(df):
  #scaler = MinMaxScaler()
  df_model = df[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'RF_distance', 'has_block', 'main_block_size']]
  df_model = df_model[(df['msa_tools'] != 'Original') & (df['RF_distance'] != -1)]
  df_model['RF_distance'] = df_model['RF_distance'] > 4
  df_model['is_AA'] = df['residue_type'] == "AA"

  df_model_features = df_model[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'has_block', 'main_block_size', 'is_AA']]
  print(df_model_features.describe())
  print(df_model_features.info())
  #print(df_model.corr().to_string())

  model = RandomForestClassifier(max_depth=3, n_estimators=100)
  X = df_model_features
  y = df_model['RF_distance']
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, train_size=0.7)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)

  print("Accuracy:", accuracy_score(y_test, predictions))


if __name__ == "__main__":
  sys.exit(main())
