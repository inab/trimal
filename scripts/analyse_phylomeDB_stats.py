
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


def main():
    df_phylomedb = pd.read_csv("/home/nicolas/Downloads/results_for_trimal.tsv", sep='\t')
    df_dessimoz = pd.read_csv("table_all_tools_AA_fixed.csv", index_col = 0)
    #print(df.head())
    df_dessimoz.loc[(df_dessimoz['msa_filter_tools'] != "None"), 'alignment_type'] = "clean"
    df_dessimoz.loc[(df_dessimoz['msa_filter_tools'] == "None"), 'alignment_type'] = "raw"

    df_phylomedb = df_phylomedb.loc[df_phylomedb['alignment_type'] != 'other', ['alignment_type', 'number of resiudes']]
    df_dessimoz = df_dessimoz.loc[:, ['alignment_type', 'num_columns']]
    df_phylomedb = df_phylomedb.rename(columns={"number of resiudes":"num_columns"})
    df_phylomedb['dataset'] = "PhylomeDB"
    df_dessimoz['dataset'] = "Tan et al."
    df = pd.concat([df_phylomedb, df_dessimoz], ignore_index=True)
    #print(df.head())

    df_phylomedb_clean_cols = df_phylomedb.loc[df_phylomedb['alignment_type'] == "clean",'num_columns']
    df_dessimoz_clean_cols = df_dessimoz.loc[df_dessimoz['alignment_type'] == "clean",'num_columns']
    print(stats.ttest_ind(a=df_phylomedb_clean_cols, b=df_dessimoz_clean_cols, equal_var=True))

    df_phylomedb_raw_cols = df_phylomedb.loc[df_phylomedb['alignment_type'] == "raw",'num_columns']
    df_dessimoz_raw_cols = df_dessimoz.loc[df_dessimoz['alignment_type'] == "raw",'num_columns']
    print(stats.ttest_ind(a=df_phylomedb_raw_cols, b=df_dessimoz_raw_cols, equal_var=True))

    '''
    q1=df["num_columns"].quantile(0.25)
    q3=df["num_columns"].quantile(0.75)
    IQR=q3-q1
    not_outliers = df.loc[(df["num_columns"] < (q3 + 1.5 * IQR)) & (df["num_columns"] > (q1 - 1.5 * IQR)), :]
    
    #df_clean_raw = df.loc[(df['alignment_type'] == "clean") | (df['alignment_type'] == "raw"), :]
    ax = sns.violinplot(x=not_outliers["alignment_type"], y=not_outliers["num_columns"], order=['raw', 'clean'], hue = not_outliers['dataset'])
    #ax = sns.boxplot(x=df["alignment_type"], y=df["num_columns"], order=['raw', 'clean'], showfliers = False, hue = df['dataset'])
    #ax = sns.swarmplot(x=not_outliers["alignment_type"], y=not_outliers["num_columns"], order=['raw', 'clean'], hue = not_outliers['dataset'], color=".25")
    ax.get_legend().set_title('')
    ax.set_xlabel("Alignment type")
    ax.set_ylabel("Number of columns")
    #sns.violinplot( x='alignment_type', y='number of resiudes', data=df)
    #plt.show()
    '''


if __name__ == "__main__":
  sys.exit(main())