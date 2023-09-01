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
from sklearn.pipeline import Pipeline
from sklearn import svm
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.inspection import DecisionBoundaryDisplay, permutation_importance
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
    parser.add_argument("--categorical_features", dest="categorical_features",
                        default=False, action="store_true", help="Use categorical features")
    parser.add_argument("--max_depth", dest="max_depth",
                        default=2, type=int, help="Tree max depth")
    parser.add_argument("--ratio", dest="ratio", default=False,
                        action="store_true", help="Add ratio parameters to the model")
    parser.add_argument("--taxon", dest="taxon", required=False, type=str,
                        choices=["Bacteria", "Eukaryotes", "Fungi"], help="Taxon")
    parser.add_argument("--blocks", dest="blocks", default=False,
                        action="store_true", help="Has some block")
    parser.add_argument("--msa_tool", dest="msa_tool", required=False, type=str,
                        choices=["ClustalW", "ClustalW2", "Mafft", "Prank", "T-Coffee"], help="MSA tool")
    parser.add_argument("--msa_filter_tool", dest="msa_filter_tool", required=False, type=str,
                        choices=["trimAl"], help="MSA tool")
    parser.add_argument("-c", "--criterion", dest="criterion", required=False, default="gini",
                        type=str, choices=["gini", "entropy"], help="The function to measure the quality of a split")
    parser.add_argument("--min_columns", dest="min_columns", required=False, type=int,
                        help="Minimin number of columns of the MSA alignment (unfiltered)")
    parser.add_argument("--min_seqs", dest="min_seqs", required=False,
                        type=int,  help="Minimin number of sequences of the alignment")
    parser.add_argument("--exclude", dest="excluded_features", required=False,
                        type=list,  help="List of features to exclude")
    parser.add_argument("--sample", dest="sample_size", required=False,
                        type=int,  help="Sample size from the dataset")
    parser.add_argument("--shapley", dest="shapley_values", default=False,
                        action="store_true",  help="Compute shapley values")
    parser.add_argument("--mdi", dest="mdi", default=False,
                        action="store_true",  help="Compute MDI values")
    parser.add_argument("--perm", dest="perm_imp", default=False,
                        action="store_true",  help="Compute permutation importance values")

    args = parser.parse_args()

    if args.explore:
        explore_data(args.max_depth, diff=args.compare, ratio=args.ratio, taxon=args.taxon,
                     tool=args.msa_tool, filter_tool=args.msa_filter_tool, criterion=args.criterion, min_columns=args.min_columns, min_seqs=args.min_seqs)
    if args.svm:
        run_svm_classifier(diff=args.compare, taxon=args.taxon,
                           tool=args.msa_tool, filter_tool=args.msa_filter_tool, min_columns=args.min_columns, min_seqs=args.min_seqs, sample_size=args.sample_size, perm_imp=args.perm_imp)
    if args.tree:
        run_decision_tree_classifier(blocks=args.blocks, use_categorical_features=args.categorical_features, max_depth=args.max_depth, diff=args.compare, taxon=args.taxon,
                                     tool=args.msa_tool, filter_tool=args.msa_filter_tool, criterion=args.criterion, min_columns=args.min_columns, min_seqs=args.min_seqs, sample_size=args.sample_size)
    if args.forest:
        run_random_forest_classifier(args.max_depth, diff=args.compare, taxon=args.taxon,
                                     tool=args.msa_tool, filter_tool=args.msa_filter_tool, criterion=args.criterion, min_columns=args.min_columns, min_seqs=args.min_seqs, sample_size=args.sample_size, shapley_values=args.shapley_values, mdi=args.mdi)
    if args.regression:
        log_regression(diff=args.compare, taxon=args.taxon,
                           tool=args.msa_tool, filter_tool=args.msa_filter_tool, min_columns=args.min_columns, min_seqs=args.min_seqs, sample_size=args.sample_size)
    if args.knn:
        run_knn_classifier(diff=args.compare, taxon=args.taxon,
                           tool=args.msa_tool, filter_tool=args.msa_filter_tool, min_columns=args.min_columns, min_seqs=args.min_seqs, sample_size=args.sample_size, shapley_values=args.shapley_values, perm_imp=args.perm_imp)
    if args.nn:
        run_nn_classifier(args.max_depth, diff=args.compare, ratio=args.ratio, taxon=args.taxon,
                           tool=args.msa_tool, filter_tool=args.msa_filter_tool, criterion=args.criterion, min_columns=args.min_columns, min_seqs=args.min_seqs)


def preprocess_data(blocks, taxon, tool, filter_tool, min_columns, min_seqs, sample_size, use_categorical_features=False):
    
    table_filename = "AA_stats.csv"
    df = pd.read_csv(table_filename, index_col=0)

    #df.loc[df["main_block_size"] < 0, "main_block_size"] = 0
    df["main_block_size"] = df["right_block_column"] - df["left_block_column"]
    df["perc_main_block_size"] = df["main_block_size"] / df["msa_columns"]
    df["right_block_column"] = df["right_block_column"] / df["msa_columns"]
    df.loc[df["right_block_column"] < 0, "right_block_column"] = -1
    df["left_block_column"] = df["left_block_column"] / df["msa_columns"]
    df.loc[df["left_block_column"] < 0, "left_block_column"] = -1
    df["blocks_diff_weighted"] = df["blocks_diff"] * df["percent_conserved_columns"]

    # clean dataset
    df = df.loc[(df['msa_tools'] != 'None') & (df['msa_filter_tools']
                                               != 'None') & ((df['RF_distance_diff'] % 2) == 0) &
                                               (df['error'] == False), :]

    df_model = df.copy()
    if taxon:
        df_model = df_model[df_model["taxon"] == taxon]
        df_model = df_model.drop(["taxon"], axis=1)
    if tool:
        df_model = df_model.loc[(df_model['msa_tools'] == tool) & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
        df_model = df_model.drop(["msa_tools"], axis=1)
    if filter_tool:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] == filter_tool) & (df_model['RF_distance'] != -1), :]
        df_model = df_model.drop(["msa_filter_tools"], axis=1)
    if tool == None and filter_tool == None:
        df_model = df_model.loc[(df_model['msa_tools'] != 'None') & (
            df_model['msa_filter_tools'] != 'None') & (df_model['RF_distance'] != -1), :]
    if min_columns:
        df_model = df_model.loc[(df_model["msa_columns"] >= min_columns), :]
    if min_seqs:
        df_model = df_model.loc[(df_model["num_sequences"] >= min_seqs), :]
    if blocks:
        df_model = df_model.loc[(df_model["num_blocks"] > 0), :]


    # Keep only computed features which are relative to the msa size and remove those highly correlated
    df_model = df_model.drop(["error", "min_columns", "max_columns", "residue_type", "RF_distance", "problem_num", "blocks_diff"], axis=1)
    df_model = df_model.drop(["main_block_size", "right_block_column",
                 "removed_columns", "gappy_columns_50"], axis=1)

    df_model.loc[(df_model['RF_distance_diff'] > 0),
                    'RF_distance_diff'] = 2
    df_model.loc[(df_model['RF_distance_diff'] == 0),
                    'RF_distance_diff'] = 1
    df_model.loc[(df_model['RF_distance_diff'] < 0),
                    'RF_distance_diff'] = 0
    # df_model.loc[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0), 'RF_distance_diff'] = 2
    # df_model = df_model.drop(df_model[(df_model['RF_distance_diff'] == 1) & (df_model['RF_distance'] == 0)].index)
 
    if use_categorical_features:
        if taxon == None:
            enc = OneHotEncoder()
            enc_df = pd.DataFrame(enc.fit_transform(
                df_model[['taxon']]).toarray())
            enc_df.columns = enc.get_feature_names_out()
            df_model = df_model.reset_index(drop=True)
            df_model = df_model.join(enc_df)
            df_model = df_model.drop(["taxon"], axis=1)

        if filter_tool == None:
            enc = OneHotEncoder()
            enc_df = pd.DataFrame(enc.fit_transform(
                df_model[['msa_filter_tools']]).toarray())
            enc_df.columns = enc.get_feature_names_out()
            df_model = df_model.reset_index(drop=True)
            df_model = df_model.join(enc_df)
            df_model = df_model.drop(["msa_filter_tools"], axis=1)

        if tool == None:
            enc = OneHotEncoder()
            enc_df = pd.DataFrame(enc.fit_transform(
                df_model[['msa_tools']]).toarray())
            enc_df.columns = enc.get_feature_names_out()
            df_model = df_model.reset_index(drop=True)
            df_model = df_model.join(enc_df)
            df_model = df_model.drop(["msa_tools"], axis=1)
    
    else:
        df_model = df_model.select_dtypes(include=['number'])

    
    class_feature = 'RF_distance_diff'
    df_model = df_model.dropna()

    if sample_size:
        df_model = resample(df_model, n_samples=sample_size, stratify=df_model[class_feature])

    print(df_model.head())
    print(df_model.describe())

    return df_model


def compute_permutation_importance(model, X_test, y_test, features):
    perm_importance = permutation_importance(model, X_test, y_test, n_jobs=4, max_samples=5, n_repeats=10)
    sorted_idx = perm_importance.importances_mean.argsort()
    plt.barh(np.array(features)[sorted_idx], perm_importance.importances_mean[sorted_idx])
    plt.xlabel("Permutation Importance")
    plt.show()


def compute_shapley_values(model, X_train, X_test, is_auto_explainer):
    #explainer = shap.Explainer(model) if is_auto_explainer else shap.KernelExplainer(model.predict_proba, shap.kmeans(X_train, 50))
    explainer = shap.explainers.Linear(model, X_train)
    shap_values = explainer.shap_values(X_test)
    shap.summary_plot(shap_values, X_test, class_names=[
                      "worse", "unchanged", "better"])
    shap.summary_plot(shap_values[0], X_test)
    shap.summary_plot(shap_values[1], X_test)
    shap.summary_plot(shap_values[2], X_test)
    shap.decision_plot(explainer.expected_value[0], shap_values[0], X_test.columns, ignore_warnings=True)
    shap.decision_plot(explainer.expected_value[1], shap_values[1], X_test.columns, ignore_warnings=True)
    shap.decision_plot(explainer.expected_value[2], shap_values[2], X_test.columns, ignore_warnings=True)


def run_decision_tree_classifier(blocks, use_categorical_features, max_depth, diff, taxon, tool, filter_tool, criterion, min_columns, min_seqs, sample_size):
    df_model = preprocess_data(blocks, taxon, tool, filter_tool, min_columns, min_seqs, sample_size, use_categorical_features)
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    features = [colname for colname in df_model.columns if colname != class_feature]
    print(features)

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
    residue_type = "AA"
    tree_filename_prex = "tree_%s_" % residue_type
    if taxon:
        tree_filename_prex += "%s_" % taxon
    if tool:
        tree_filename_prex += "%s_" % tool
    if diff:
        tree_filename_prex += "diff_"
    if min_columns:
        tree_filename_prex += "%s_min_columns_" % min_columns
    if min_seqs:
        tree_filename_prex += "%s_min_seqs_" % min_seqs
    tree_filename = tree_filename_prex + \
        str(max_depth) + "_" + criterion + ".png"
    print("saved " + tree_filename)
    graph.write_png(tree_filename)


def run_random_forest_classifier(max_depth, diff, taxon, tool, filter_tool, criterion, min_columns, min_seqs, sample_size, shapley_values, mdi):
    df_model = preprocess_data(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size)
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    features = [colname for colname in df_model.columns if colname != class_feature]
    print(features)

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

    ConfusionMatrixDisplay.from_estimator(model, X_test, y_test)
    plt.show()

    if shapley_values: compute_shapley_values(model, X_train, X_test, is_auto_explainer=True)
    if mdi:
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


def run_svm_classifier(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size, perm_imp):
    df_model = preprocess_data(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size)
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    features = [colname for colname in df_model.columns if colname != class_feature]
    print(features)

    X = df_model[features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8)

    print(len(X_train))
    print(len(X_test))

    model = Pipeline([('scale', StandardScaler()), ('svm', svm.SVC(
        decision_function_shape='ovr', class_weight='balanced', cache_size=3500, verbose=True, break_ties=True))])
    model.fit(X_train, y_train)

    predictions = model.predict(X_test)
    print(classification_report(y_test, predictions))

    ConfusionMatrixDisplay.from_estimator(model, X_test, y_test)
    plt.show()

    if perm_imp: compute_permutation_importance(model, X_test, y_test, features)
        

def run_knn_classifier(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size, shapley_values, perm_imp):
    df_model = preprocess_data(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size)
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    features = [colname for colname in df_model.columns if colname != class_feature]
    print(features)

    X = df_model[features]
    y = df_model[class_feature]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8, stratify=y)

    print(len(X_train))
    print(len(X_test))

    
    model = Pipeline([('scale', StandardScaler()), ('knn', KNeighborsClassifier(n_neighbors=3, n_jobs=4))])
    model.fit(X_train, y_train)

    predictions = model.predict(X_test)
    print(classification_report(y_test, predictions))

    ConfusionMatrixDisplay.from_estimator(model, X_test, y_test)
    plt.show()

    # if shapley_values: compute_shapley_values(model['knn'], X_train, X_test, is_auto_explainer=False)
    if perm_imp: compute_permutation_importance(model, X_test, y_test, features)


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
    # model = make_pipeline(StandardScaler(), GridSearchCV(mlp, parameter_space, n_jobs=4, cv=3, verbose=True))

    # model = GridSearchCV(mlp, parameter_space, n_jobs=4, cv=3, verbose=True)
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


def log_regression(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size):
    df_model = preprocess_data(diff, taxon, tool, filter_tool, min_columns, min_seqs, sample_size)
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    features = [colname for colname in df_model.columns if colname != class_feature]
    print(features)

    model = LogisticRegression(class_weight='balanced')
    X = df_model[features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, shuffle=True, train_size=0.8)
    
    std_scaler = StandardScaler()
    X_train_scaled = pd.DataFrame(std_scaler.fit_transform(X_train), columns=X_train.columns)
    X_test_scaled = pd.DataFrame(std_scaler.transform(X_test), columns=X_test.columns)
    

    model.fit(X_train_scaled, y_train)
    print(len(X_train_scaled))
    print(len(X_test_scaled))

    cdf = pd.concat([pd.DataFrame(X.columns), pd.DataFrame(
        np.transpose(model.coef_))], axis=1)
    print(cdf)

    predictions = model.predict(X_test_scaled)

    mae = mean_absolute_error(y_test, predictions)
    mse = mean_squared_error(y_test, predictions)
    r2 = r2_score(y_test, predictions)

    print(classification_report(y_test, predictions))

    print("The model performance for testing set")
    print("--------------------------------------")
    print('MAE is {}'.format(mae))
    print('MSE is {}'.format(mse))
    print('R2 score is {}'.format(r2))

    compute_shapley_values(model, X_train_scaled, X_test_scaled, is_auto_explainer=True)


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
    "avg_seq_identity_diff_weighted", "RF_distance_diff"]
    class_feature = 'RF_distance_diff' if diff else 'RF_distance'
    class_feature = "taxon"
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
    

    df_model = resample(df_model, n_samples=100000,
        stratify=df_model[class_feature])
    
    X = df_model[numerical_features]
    y = df_model[class_feature]
    X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle = True, train_size = 0.8)
    
    sc = StandardScaler()
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)


    import plotly.express as px

    from MulticoreTSNE import MulticoreTSNE as TSNE


    print(df_model["taxon"].head())


    #tsne = TSNE(n_components=2, random_state=42, n_jobs=4)
    #X_tsne = tsne.fit_transform(X)
    
    tsne = TSNE(n_components=2, perplexity=5, n_jobs=4, verbose=1)
    X_tsne = tsne.fit_transform(X)
    print("tsne computed")
    print(tsne.kl_divergence_)

    fig = px.scatter(x=X_tsne[:, 0], y=X_tsne[:, 1], color=df_model["taxon"])
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
