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

import os
import sys
import argparse
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

def main():
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


if __name__ == "__main__":
  sys.exit(main())
