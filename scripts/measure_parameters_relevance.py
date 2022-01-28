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

from ctypes import alignment
import os
import sys
import argparse
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import MinMaxScaler
import random

def main():
  '''
  df = pd.read_csv('https://raw.githubusercontent.com/datagy/data/main/insurance.csv')
  print(df.head())
  print(df.info())
  print(df.corr())

  model = LinearRegression()
  non_smokers = df[df['smoker'] == 'no']
  X = non_smokers[['age', 'bmi']]
  y = non_smokers['charges']
  X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=True, train_size=0.3)
  model.fit(X_train, y_train)

  predictions = model.predict(X_test)
  r2 = r2_score(y_test, predictions)
  rmse = mean_squared_error(y_test, predictions, squared=False)

  print('The r2 is: ', r2)
  print('The rmse is: ', rmse)
  print(model.coef_)
  print(model.intercept_)
  '''

  '''
  alignment_index = 0
  with open(r'test_working_files/gap_stats.txt', 'r') as infile, \
     open(r'test_working_files/gap_stats_clean.txt', 'w') as outfile:
    data = infile.read()
    data = data.replace("\t\t", "\t")
    outfile.write(data)

  df = pd.read_table("test_working_files/gap_stats_clean.txt", index_col=False, skiprows=3, names=["#_residues",
    "%_length", "cumulative_# residues", "%_cumulative_length", "gaps/column", "%_gaps/column", "gap_score/column"])
  df["alingment_index"] = alignment_index
  
  print(df.head())
  print(df.info())

  num_colums = sum(df["#_residues"])
  '''

  
  num_sequences = pd.read_table("test_working_files/number_sequences.txt", header=None, names=["num_sequences"])
  num_blocks = pd.read_table("test_working_files/blocks_outputs.txt", header=None, names=["num_blocks"])
  num_columns = pd.read_table("test_working_files/number_columns.txt", header=None, names=["num_columns"])
  avg_seq_identity = pd.read_table("test_working_files/avg_seq_identity.txt", header=None, names=["avg_seq_identity"])
  avg_seq_overlap = pd.read_table("test_working_files/avg_seq_overlap.txt", header=None, names=["avg_seq_overlap"])
  num_identical_columns = pd.read_table("test_working_files/number_identical_columns.txt", header=None, names=["num_identical_columns"])

  df = pd.concat([num_sequences, num_blocks, num_columns, avg_seq_identity, avg_seq_overlap,
    num_identical_columns], axis=1)
  df["score"] = [random.random() for i in range(df.shape[0])]


  print(df.head())
  print(df.info())
  print(df.corr().to_string())

  #run_model(df)


def run_model(df):
  scaler = MinMaxScaler()
  df = pd.DataFrame(scaler.fit_transform(df),
            columns=df.columns, index=df.index) 

  model = LinearRegression()
  X = df[['num_sequences', 'num_blocks', 'num_columns', 'avg_seq_identity', 'avg_seq_overlap', 'num_identical_columns']]
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
