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

import sys
import pandas as pd
from joblib import Parallel, delayed


def main():

  df = pd.read_csv("RF_distances_AA_automated3.csv")
  df = df.loc[(df["MSA_filter_tool"] == "trimAl_gappyout") | (df["MSA_filter_tool"] == "trimAl_automated2_gaps_30_100") | (df["MSA_filter_tool"] == "trimAl_automated2_gaps_30_250"), :]
  add_automated3(df, 123)


def add_automated3(df, problem_num):
  for taxon in ["Bacteria", "Eukaryotes", "Fungi"]:
    df_problem = df.loc[(df["problem_num"] == problem_num) & (df["taxon"] == taxon), :]
    if df_problem.loc[(df_problem["MSA_filter_tool"] == "trimAl_gappyout"), "number_columns"].values[0] >\
       df_problem.loc[(df_problem["MSA_filter_tool"] == "trimAl_automated2_gaps_30_100"), "number_columns"].values[0]:
        df_problem = df_problem.loc[(df_problem["MSA_filter_tool"] == "trimAl_gappyout"), :]
    else:
      df_problem = df_problem.loc[(df_problem["MSA_filter_tool"] == "trimAl_gappyout"), :]
    
    df_problem["MSA_filter_tool"] = "trimAl_automated3_1"
    print(df_problem.to_string())




if __name__ == "__main__":
  sys.exit(main())
