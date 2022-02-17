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
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import MinMaxScaler

def main():
  num_sequences = pd.read_table("test_working_files/number_sequences.txt", header=None, names=["num_sequences"], dtype="int")
  num_blocks = pd.read_table("test_working_files/blocks_outputs.txt", header=None, names=["num_blocks"], dtype="int")
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

  df = pd.concat([num_sequences, num_blocks, num_columns, min_columns, max_columns, avg_seq_identity, avg_gaps, RF_distance, residue_type, taxon, problem_num, error, msa_tools, msa_filter_tools], axis=1)
  df["has_block"] = df["num_blocks"] > 0

  print(df.head())
  print(df.tail())
  print(df.info())

  with open('table.csv', 'w') as file:
    file.write(df.to_csv())

  '''
  print(df.corr().to_string())

  with open('table.html', 'w') as file:
    file.write(df.to_html())

  with open('table.csv', 'w') as file:
    file.write(df.to_csv())
  '''

  #run_model(df)


def run_model(df):
  scaler = MinMaxScaler()
  df = pd.DataFrame(scaler.fit_transform(df),
            columns=df.columns, index=df.index) 

  model = LinearRegression()
  X = df[['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps']]
  y = df['score']
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, train_size=0.3)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)
  r2 = r2_score(y_test, predictions)
  rmse = mean_squared_error(y_test, predictions, squared=False)

  print('The r2 is: ', r2)
  print('The rmse is: ', rmse)
  print(model.coef_)
  print(model.intercept_)


if __name__ == "__main__":
  sys.exit(main())
