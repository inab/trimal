#!/usr/bin/python

#
#  'measure_parameters_relevance.py'
#
#   Script implemented to work with trimAl to analyze parameters statistics
#   and decide which are more relevant to the quality of the alignment.
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
import os
import pandas as pd
import pydotplus
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import shap
import pingouin as pg
import joblib

from sklearn.model_selection import cross_validate, train_test_split
from sklearn.metrics import accuracy_score, mean_absolute_error, mean_squared_error, r2_score, classification_report, ConfusionMatrixDisplay
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import resample
from sklearn.pipeline import make_pipeline
from sklearn import svm
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV

from matplotlib.colors import ListedColormap


df = pd.DataFrame()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--tree", dest="tree", default=False,
                        action="store_true", help="Build decision tree")
    parser.add_argument("--forest", dest="forest", default=False,
                        action="store_true", help="Build random forest")
    parser.add_argument("--regression", dest="regression", default=False,
                        action="store_true", help="Build linear regression")
    parser.add_argument("--svm", dest="svm", default=False,
                        action="store_true", help="Build svm")
    parser.add_argument("--knn", dest="knn", default=False,
                        action="store_true", help="Build knn")
    parser.add_argument("--nn", dest="nn", default=False,
                        action="store_true", help="Build nn")
    parser.add_argument("--explore", dest="explore", default=False,
                        action="store_true", help="Explore dataset")
    parser.add_argument("--compare", dest="compare", default=False,
                        action="store_true", help="Compare filtered results with original MSA")
    parser.add_argument("--recalculate", dest="recalculate",
                        default=False, action="store_true", help="Recalculate values")
    parser.add_argument("--max_depth", dest="maxDepth",
                        default=2, type=int, help="Tree max depth")
    parser.add_argument("--ratio", dest="ratio", default=False,
                        action="store_true", help="Add ratio parameters to the model")
    parser.add_argument("-t", "--type", dest="residueType", required=True,
                        type=str, choices=["AA", "DNA"], help="Residue type")
    parser.add_argument("--taxon", dest="taxon", required=False, type=str,
                        choices=["Bacteria", "Eukaryotes", "Fungi"], help="Taxon")
    parser.add_argument("--msa_tool", dest="msaTool", required=False, type=str,
                        choices=["ClustalW", "ClustalW2", "Mafft", "Prank", "T-Coffee"], help="MSA tool")
    parser.add_argument("-c", "--criterion", dest="criterion", required=False, default="gini",
                        type=str, choices=["gini", "entropy"], help="The function to measure the quality of a split")
    parser.add_argument("--min_columns", dest="minColumns", required=False, type=int,
                        help="Minimin number of columns of the MSA alignment (unfiltered)")
    parser.add_argument("--min_seqs", dest="minSeqs", required=False,
                        type=int,  help="Minimin number of sequences of the alignment")

    args = parser.parse_args()

    global df
    table_filename = "AA_stats.csv"
    df = pd.read_csv(table_filename, index_col=0)

    if args.recalculate and args.ratio:
        df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
        df.loc[df["main_block_size"] < 0, "main_block_size"] = 0
        df["columns/sequence"] = df["msa_columns"] / df["num_sequences"]
        df["blocks/columns"] = df["num_blocks"] / df["msa_columns"]
        df["perc_main_block_size"] = df["main_block_size"] / df["msa_columns"]

    df["has_block"] = df["num_blocks"] > 0
    df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
    df["perc_main_block_size"] = df["main_block_size"] / df["msa_columns"]
    df["right_block_column"] = df["right_block_column"] / df["msa_columns"]
    df.loc[df["right_block_column"] < 0, "right_block_column"] = -1
    df["left_block_column"] = df["left_block_column"] / df["msa_columns"]
    df.loc[df["left_block_column"] < 0, "left_block_column"] = -1 # consider only msas with some block to see the influence of this variable (?)

    # clean dataset
    df = df.loc[(df['msa_tools'] != 'None') & (df['msa_filter_tools']
                                               != 'None') & ((df['RF_distance_diff'] % 2) == 0) &
                                               (df['error'] == False), :]
    df = df.drop(["error", "min_columns", "max_columns"], axis=1)

    # Keep only computed features which are relative to the msa size and remove those highly correlated
    df = df.drop(["main_block_size", "right_block_column",
                 "removed_columns", "gappy_columns_50"], axis=1)

    # Impute some data or discard (e.g. left_block_column when it's -1)?

    # Repeat with four classes?

    if args.explore:
        explore_data(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                     tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)
    if args.svm:
        run_svm_classifier(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                           tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)
    if args.tree:
        run_decision_tree_classifier(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                                     tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)
    if args.forest:
        run_random_forest_classifier(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                                     tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)
    if args.regression:
        log_regression(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                       tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)
    if args.knn:
        run_knn_classifier(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                           tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)
    if args.nn:
        run_nn_classifier(df, args.maxDepth, diff=args.compare, ratio=args.ratio, residue_type=args.residueType, taxon=args.taxon,
                           tool=args.msaTool, criterion=args.criterion, min_columns=args.minColumns, min_seqs=args.minSeqs)


def run_decision_tree_classifier(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
    features = ['num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted',
                "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += ['percent_conserved_columns']
        # Ignore removed columns because is highly correlated with msa_columns and depends on the size of the msa
    if ratio:
        features += ['columns/sequence', 'blocks/columns']
    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    else:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if diff:
        df_model.loc[(df_model['RF_distance_diff'] > 0),
                     'RF_distance_diff'] = 2
        df_model.loc[(df_model['RF_distance_diff'] == 0),
                     'RF_distance_diff'] = 1
        df_model.loc[(df_model['RF_distance_diff'] < 0),
                     'RF_distance_diff'] = 0
        # df_model = df_model.drop(df_model[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0)].index)
    else:
        df_model['RF_distance'] = df_model['RF_distance'] < 4

    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(
        df_model[['msa_filter_tools']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()

    if not tool:
        enc = OneHotEncoder()
        enc_df = pd.DataFrame(enc.fit_transform(
            df_model[['msa_tools']]).toarray())
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

    # df_model = resample(df_model, n_samples=70000, stratify=df_model[class_feature])

    X = df_model[features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8)

    model = DecisionTreeClassifier(
        max_depth=max_depth, criterion=criterion, class_weight='balanced')
    model.fit(X_train, y_train)

    predictions = model.predict(X_test)

    print(classification_report(y_test, predictions))

    # class_names = ['worse', 'unchanged-wrong', 'better', 'unchanged-right'] if diff else ['different', 'similar']
    class_names = ['worse', 'unchanged', 'better'] if diff else [
        'different', 'similar']

    dot_data = export_graphviz(model, filled=True, rounded=True, special_characters=True, proportion=True, precision=2,
                               feature_names=df_model[features].columns, class_names=class_names)
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
    tree_filename = tree_filename_prex + \
        str(max_depth) + "_" + criterion + ".png"
    print("saved " + tree_filename)
    graph.write_png(tree_filename)


def run_random_forest_classifier(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
    features = ['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted',
                "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += ['percent_conserved_columns']
    if ratio:
        features += ['columns/sequence', 'blocks/columns']
    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    else:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if diff:
        df_model.loc[(df_model['RF_distance_diff'] > 0),
                     'RF_distance_diff'] = 2
        df_model.loc[(df_model['RF_distance_diff'] == 0),
                     'RF_distance_diff'] = 1
        df_model.loc[(df_model['RF_distance_diff'] < 0),
                     'RF_distance_diff'] = 0
        # df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 2
        # df_model = df_model.drop(df_model[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0)].index)
    else:
        df_model['RF_distance'] = df_model['RF_distance'] < 4

    # add parameter to ignore
    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(
        df_model[['msa_filter_tools']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()

    # add parameter to ignore
    if not tool:
        enc = OneHotEncoder()
        enc_df = pd.DataFrame(enc.fit_transform(
            df_model[['msa_tools']]).toarray())
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

    # df_model = resample(df_model, n_samples=70000, stratify=df_model[class_feature])

    X = df_model[features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8)

    print(len(X_train))
    print(len(X_test))

    model = RandomForestClassifier(
        max_depth=max_depth, criterion=criterion, n_jobs=4, min_samples_leaf=10, class_weight='balanced')
    model.fit(X_train, y_train)

    predictions = model.predict(X_test)

    print(classification_report(y_test, predictions))

    return
    explainer = shap.Explainer(model)
    shap_values = explainer.shap_values(X_test)
    shap.summary_plot(shap_values, X_test, class_names=[
                      "worse", "unchanged", "better"])
    # shap.dependence_plot("num_columns", shap_values[0], X_test,interaction_index="percent_conserved_columns")

    ConfusionMatrixDisplay.from_estimator(model, X_test, y_test)
    plt.show()

    importances = model.feature_importances_
    std = np.std(
        [tree.feature_importances_ for tree in model.estimators_], axis=0)

    forest_importances = pd.Series(importances, index=features)

    fig, ax = plt.subplots()
    forest_importances.plot.bar(yerr=std, ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()
    plt.show()


def run_svm_classifier(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
    features = ['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted',
                "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += ['percent_conserved_columns']
    if ratio:
        features += ['columns/sequence', 'blocks/columns']
    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    else:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if diff:
        df_model.loc[(df_model['RF_distance_diff'] > 0),
                     'RF_distance_diff'] = 2
        df_model.loc[(df_model['RF_distance_diff'] == 0),
                     'RF_distance_diff'] = 1
        df_model.loc[(df_model['RF_distance_diff'] < 0),
                     'RF_distance_diff'] = 0
        # df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 2
    else:
        df_model['RF_distance'] = df_model['RF_distance'] < 4

    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(
        df_model[['msa_filter_tools']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()

    if not tool:
        enc = OneHotEncoder()
        enc_df = pd.DataFrame(enc.fit_transform(
            df_model[['msa_tools']]).toarray())
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
    df_model = resample(df_model, n_samples=70000,
                        stratify=df_model[class_feature])

    X = df_model[features]
    y = df_model[class_feature]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8)

    print(len(X_train))
    print(len(X_test))

    '''
  features = [ 'num_columns', 'percent_conserved_columns']
  _, ax = plt.subplots()
  cmap_light = ListedColormap(["orange", "cyan", "cornflowerblue"])
  cmap_bold = ["darkorange", "c", "darkblue"]
  DecisionBoundaryDisplay.from_estimator(
        model,
        X_train,
        cmap=cmap_light,
        ax=ax,
        response_method="predict",
        plot_method="pcolormesh",
        xlabel="num_columns",
        ylabel="percent_conserved_columns",
        shading="auto",
    )

  sns.scatterplot(
        x=X_train["num_columns"],
        y=X_train["percent_conserved_columns"],
        hue=y_train,
        palette=cmap_bold,
        alpha=1.0,
        edgecolor="black",
    )

  plt.show()
  '''

    model = make_pipeline(StandardScaler(), svm.SVC(
        decision_function_shape='ovr', class_weight='balanced', cache_size=3500, verbose=True))
    # model = svm.SVC(decision_function_shape='ovr', class_weight='balanced', cache_size=3500, verbose=True)
    model.fit(X_train, y_train)

    predictions = model.predict(X_test)

    print(classification_report(y_test, predictions))


def run_knn_classifier(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
    features = ['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted',
                "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += ['percent_conserved_columns']
    if ratio:
        features += ['columns/sequence', 'blocks/columns']
    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    else:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if diff:
        df_model.loc[(df_model['RF_distance_diff'] > 0),
                     'RF_distance_diff'] = 2
        df_model.loc[(df_model['RF_distance_diff'] == 0),
                     'RF_distance_diff'] = 1
        df_model.loc[(df_model['RF_distance_diff'] < 0),
                     'RF_distance_diff'] = 0
        # df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 2
        # df_model = df_model.drop(df_model[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0)].index)
    else:
        df_model['RF_distance'] = df_model['RF_distance'] < 4

    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(
        df_model[['msa_filter_tools']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()

    if not tool:
        enc = OneHotEncoder()
        enc_df = pd.DataFrame(enc.fit_transform(
            df_model[['msa_tools']]).toarray())
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
    # df_model = resample(df_model, n_samples=30000,
     #                   stratify=df_model[class_feature])

    X = df_model[features]
    y = df_model[class_feature]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8, stratify=y)

    print(len(X_train))
    print(len(X_test))

    '''
  features = [ 'num_columns', 'percent_conserved_columns']
  _, ax = plt.subplots()
  cmap_light = ListedColormap(["orange", "cyan", "cornflowerblue"])
  cmap_bold = ["darkorange", "c", "darkblue"]
  DecisionBoundaryDisplay.from_estimator(
        model,
        X_train,
        cmap=cmap_light,
        ax=ax,
        response_method="predict",
        plot_method="pcolormesh",
        xlabel="num_columns",
        ylabel="percent_conserved_columns",
        shading="auto",
    )

  sns.scatterplot(
        x=X_train["num_columns"],
        y=X_train["percent_conserved_columns"],
        hue=y_train,
        palette=cmap_bold,
        alpha=1.0,
        edgecolor="black",
    )

  plt.show()
  '''

    model = make_pipeline(
        StandardScaler(), KNeighborsClassifier(n_neighbors=3, n_jobs=-1))  # scales both train and test data
    model.fit(X_train, y_train)

    predictions = model.predict(X_test)

    print(classification_report(y_test, predictions))


def run_nn_classifier(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
    features = ['num_sequences', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted',
                "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += ['percent_conserved_columns']
    if ratio:
        features += ['columns/sequence', 'blocks/columns']
    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    else:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if diff:
        df_model.loc[(df_model['RF_distance_diff'] > 0),
                     'RF_distance_diff'] = 2
        df_model.loc[(df_model['RF_distance_diff'] == 0),
                     'RF_distance_diff'] = 1
        df_model.loc[(df_model['RF_distance_diff'] < 0),
                     'RF_distance_diff'] = 0
        # df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 2
    else:
        df_model['RF_distance'] = df_model['RF_distance'] < 4

    enc = OneHotEncoder()
    enc_df = pd.DataFrame(enc.fit_transform(
        df_model[['msa_filter_tools']]).toarray())
    enc_df.columns = enc.get_feature_names_out()
    df_model = df_model.reset_index(drop=True)
    df_model = df_model.join(enc_df)
    features += enc.get_feature_names_out().tolist()

    if not tool:
        enc = OneHotEncoder()
        enc_df = pd.DataFrame(enc.fit_transform(
            df_model[['msa_tools']]).toarray())
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
    # df_model = resample(df_model, n_samples=30000,
     #                   stratify=df_model[class_feature])

    X = df_model[features]
    y = df_model[class_feature]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8, stratify=y)

    print(len(X_train))
    print(len(X_test))

    '''
    parameter_space = {
        'hidden_layer_sizes': [(50, 20, 10), (10, 20, 10), (100, 20, 10), (30, 50, 25)],
        'activation': ['tanh', 'relu'],
        'solver': ['sgd', 'adam'],
        'alpha': [0.0001, 0.05],
        'learning_rate': ['constant', 'adaptive'],
    }
    '''

    model = MLPClassifier(verbose=True, max_iter=100, activation='relu', alpha=0.05,
                          hidden_layer_sizes=(100, 20, 10), learning_rate='constant', solver='adam')

    # Scale the features using StandardScaler
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    # model = make_pipeline(StandardScaler(), GridSearchCV(mlp, parameter_space, n_jobs=-1, cv=3, verbose=True))

    # model = GridSearchCV(mlp, parameter_space, n_jobs=-1, cv=3, verbose=True)
    model.fit(X_train_scaled, y_train)

    predictions = model.predict(X_test_scaled)

    print(classification_report(y_test, predictions))

    return

    # Best paramete set
    print('Best parameters found:\n', model.best_params_)

    # All results
    means = model.cv_results_['mean_test_score']
    stds = model.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, model.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))

    predictions = model.predict(X_test_scaled)

    print(classification_report(y_test, predictions))


def log_regression(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):
    features = ['num_sequences', 'msa_columns', 'num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'perc_main_block_size', 'avg_gaps_diff_weighted',
                "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += ['percent_conserved_columns']
    if ratio:
        features += ['columns/sequence', 'blocks/columns']
    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    else:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if diff:
        df_model.loc[(df_model['RF_distance_diff'] > 0),
                     'RF_distance_diff'] = 2
        df_model.loc[(df_model['RF_distance_diff'] == 0),
                     'RF_distance_diff'] = 1
        df_model.loc[(df_model['RF_distance_diff'] < 0),
                     'RF_distance_diff'] = 0
    else:
        df_model['RF_distance'] = df_model['RF_distance'] < 4

    '''
  enc = OneHotEncoder()
  enc_df = pd.DataFrame(enc.fit_transform(
      df_model[['MSA_filter_tool']]).toarray())
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
    df_model.loc[:, features] = std_scaler.fit_transform(
        df_model.loc[:, features])

    print(df_model[features].head())
    print(df_model[features].describe())
    print(df_model[features].info())

    print(df_model[class_feature].info())

    model = LogisticRegression(class_weight='balanced')
    X = df_model[features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8)
    model.fit(X_train, y_train)
    print(len(X_train))
    print(len(X_test))

    print(model.score(X_train, y_train))

    cdf = pd.concat([pd.DataFrame(X.columns), pd.DataFrame(
        np.transpose(model.coef_))], axis=1)
    print(cdf)

    predictions = model.predict(X_test)

    mae = mean_absolute_error(y_test, predictions)
    mse = mean_squared_error(y_test, predictions)
    r2 = r2_score(y_test, predictions)

    print(classification_report(y_test, predictions))

    print("The model performance for testing set")
    print("--------------------------------------")
    print('MAE is {}'.format(mae))
    print('MSE is {}'.format(mse))
    print('R2 score is {}'.format(r2))


def explore_data(df, max_depth, diff, ratio, residue_type, taxon, tool, criterion, min_columns, min_seqs):

# Correlation matrix
    '''
    print(df.head())
    print(df.describe())
    print(df.info())

    numerical_data = df.select_dtypes(include=['number'])
    print(numerical_data.info())
    corr_matrix = numerical_data.corr().round(2)
    threshold = 0.85
    highly_correlated_features = np.where(np.abs(corr_matrix) > threshold)
    for feature_a, feature_b in zip(*highly_correlated_features):
      if feature_a != feature_b and feature_a < feature_b:
          print(
              f"Features '{numerical_data.columns[feature_a]}' and '{numerical_data.columns[feature_b]}' are highly correlated: {corr_matrix.iloc[feature_a, feature_b]}")



    sns.heatmap(corr_matrix, annot=True, vmax=1,
                vmin=-1, center=0, cmap='vlag')
    plt.show()
    '''

    # RF change plots
    '''
    colors = {
        "worse": "darkred",
        "unchanged": "darkorange",
        "better": "darkgreen"
    }

    df.loc[:, "RF_change"] = "None"
    df.loc[(df['RF_distance_diff'] > 0), 'RF_change'] = "better"
    df.loc[(df['RF_distance_diff'] == 0), 'RF_change'] = "unchanged"
    df.loc[(df['RF_distance_diff'] < 0), 'RF_change'] = "worse"

    sns.histplot(data=df, x="percent_conserved_columns", hue="RF_change")
    plt.show()

    sns.violinplot(data=df, x="RF_change",
                   y="percent_conserved_columns", palette=colors)
    plt.xlabel("RF change")
    plt.ylabel("% of conserved columns")
    plt.show()

    sns.boxplot(data=df, x="RF_change",
                y="percent_conserved_columns", palette=colors)
    plt.xlabel("RF change")
    plt.ylabel("% of conserved columns")
    plt.show()

    df_1000 = df.loc[df['num_columns'] < 1000,]

    sns.violinplot(data=df_1000, x="RF_change",
                   y="num_columns", palette=colors)
    plt.xlabel("RF change")
    plt.ylabel("Columns")
    plt.show()

    sns.boxplot(data=df_1000, x="RF_change", y="num_columns", palette=colors)
    plt.xlabel("RF change")
    plt.ylabel("Columns")
    plt.show()
    '''

    '''
    balance data
    df_int = df_model[df_model['RF_distance_diff'] == 0]
    df_majority = df_model[df_model['RF_distance_diff'] == 1]
    df_minority = df_model[df_model['RF_distance_diff'] == 2]

    # Upsample minority class
    df_minority_upsampled = resample(df_minority,
                                    replace=True,     # sample with replacement
                                    # to match majority class
                                    n_samples=len(df_int),
                                    random_state=1234) # reproducible results

    # Combine majority class with upsampled minority class
    df_upsampled = pd.concat([df_int, df_majority, df_minority_upsampled])

    # Display new class counts
    df_upsampled.value_counts()
    df_model = df_upsampled.copy()

    '''

    '''
  features = ['num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'main_block_size', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted', \
    "avg_seq_identity_diff_weighted"]
  class_feature = 'RF_distance_diff' if diff else 'RF_distance'
  if diff:
    features += [ 'percent_conserved_columns']
    # Ignore removed columns because is highly correlated with msa_columns and depends on the size of the msa
  if ratio:
    features += ['columns/sequence', 'blocks/columns']
  df_model = df.copy()
  if taxon:
    df_model = df_model[df_model["taxon"] == taxon]
  if tool:
    df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
        df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
  else:
    df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
        df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
  if min_columns:
    df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
  if min_seqs:
    df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
  if diff:
    df_model.loc[(df_model['RF_distance_diff'] > 0), 'RF_distance_diff'] = 2
    df_model.loc[(df_model['RF_distance_diff'] == 0), 'RF_distance_diff'] = 1
    df_model.loc[(df_model['RF_distance_diff'] < 0), 'RF_distance_diff'] = 0
    #df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 3
  else:
    df_model['RF_distance'] = df_model['RF_distance'] < 4

  enc = OneHotEncoder()
  enc_df = pd.DataFrame(enc.fit_transform(
      df_model[['msa_filter_tools']]).toarray())
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
  print(str(df_model))

  numerical_features = ['num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'main_block_size', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted', \
    "avg_seq_identity_diff_weighted"]

  X = df_model[numerical_features]
  y = df_model[class_feature]
  X_train, X_test, y_train, y_test = train_test_split(
      X, y, shuffle = True, train_size = 0.8)

  sc = StandardScaler()
  X_train = sc.fit_transform(X_train)
  X_test = sc.transform(X_test)

  pca = PCA(n_components=4)
  X_train = pca.fit_transform(X_train)
  X_test = pca.transform(X_test)

  print(pca.explained_variance_ratio_)
  print(abs( pca.components_ ))

  plt.xlim(-1,1)
  plt.ylim(-1,1)
  plt.xlabel("PC{}".format(1))
  plt.ylabel("PC{}".format(2))
  plt.grid()

  score = X_train[:,0:2]
  coeff = np.transpose(pca.components_[0:2, :])
  labels=None
  xs = score[:,0]
  ys = score[:,1]
  n = coeff.shape[0]
  scalex = 1.0/(xs.max() - xs.min())
  scaley = 1.0/(ys.max() - ys.min())
  plt.scatter(xs * scalex,ys * scaley, c = y_train)
  for i in range(n):
      plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
      if labels is None:
          plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+ \
                   str(i+1), color = 'g', ha = 'center', va = 'center')
      else:
          plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i],
                   color = 'g', ha = 'center', va = 'center')


  plt.show()

  var = np.cumsum(np.round(pca.explained_variance_ratio_, decimals=3)*100)
  print(var)
  plt.ylabel('% Variance Explained')
  plt.xlabel('# of Features')
  plt.title('PCA Analysis')
  plt.ylim(30,100.5)
  plt.style.context('seaborn-whitegrid')
  plt.plot(var)
  plt.show()

  import plotly.express as px

  fig = px.scatter(x=X_train[:, 0], y=X_train[:, 1], color=y_train)
  fig.update_layout(
      title="PCA visualization of Custom Classification dataset",
      xaxis_title="First Principal Component",
      yaxis_title="Second Principal Component",
  )
  fig.show()

  plt.figure(figsize=(16,10))
  sns.scatterplot(
    x=X_train[:, 0], y=X_train[:, 1],
    hue=y_train,
    palette=sns.color_palette("hls", 3),
    legend="full"
  )

  plt.show()

  from mpl_toolkits.mplot3d import Axes3D

  fig = plt.figure(figsize=(9,9))
  axes = Axes3D(fig)
  axes.set_title('PCA Representation', size=14)
  axes.set_xlabel('PC1')
  axes.set_ylabel('PC2')
  axes.set_zlabel('PC3')

  axes.scatter(X_train[:, 0], X_train[:, 1], X_train[:, 2],
               c=y_train, cmap = 'prism', s=10)
  plt.show()
  '''

    features = ['num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'has_block', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted',
    "avg_seq_identity_diff_weighted"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    if diff:
        features += [ 'percent_conserved_columns']
        # Ignore removed columns because is highly correlated with msa_columns and depends on the size of the msa
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
        #df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 3
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
    print(str(df_model))

    numerical_features = ['num_blocks', 'num_columns', 'avg_gaps', 'avg_seq_identity', 'msa_columns', 'perc_main_block_size', 'avg_gaps_diff_weighted', \
        "avg_seq_identity_diff_weighted"]
    
    X = df_model[numerical_features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle = True, train_size = 0.8)
    
    sc = StandardScaler()
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)


    import plotly.express as px

    from MulticoreTSNE import MulticoreTSNE as TSNE





    #tsne = TSNE(n_components=2, random_state=42, n_jobs=-1)
    #X_tsne = tsne.fit_transform(X)
    
    tsne = TSNE(n_components=2, perplexity=12, n_jobs=8)
    X_tsne = tsne.fit_transform(X)
    print("tsne computed")
    print(tsne.kl_divergence_)

    fig = px.scatter(x=X_tsne[:, 0], y=X_tsne[:, 1], color=y)
    fig.update_layout(
        title="t-SNE visualization of Custom Classification dataset",
        xaxis_title="First t-SNE",
        yaxis_title="Second t-SNE",
    )
    fig.show()

    return

'''
  perplexity = np.arange(115, 135, 5)
  divergence = []

  for i in perplexity:
      model = TSNE(n_components=2, perplexity=i, n_jobs=8)
      X_tsne = model.fit_transform(X_train)
      divergence.append(model.kl_divergence_)
  fig = px.line(x=perplexity, y=divergence, markers=True)
  fig.update_layout(xaxis_title="Perplexity Values", yaxis_title="Divergence")
  fig.update_traces(line_color="red", line_width=1)
  fig.show()

  

  import umap
  
  reducer = umap.UMAP()
  embedding = reducer.fit_transform(scaled_X_train)

  plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=[sns.color_palette()[x] for x in y_train])
  plt.gca().set_aspect('equal', 'datalim')
  plt.title('UMAP projection', fontsize=24)
  plt.show()



  explainer = shap.Explainer(model)
  shap_values = explainer.shap_values(X_test)
  shap.summary_plot(shap_values, X_test)
  shap.dependence_plot("num_columns", shap_values[0], X_test,interaction_index="percent_conserved_columns")
  
    # from sklearn.manifold import TSNE
    '''


if __name__ == "__main__":
    sys.exit(main())
