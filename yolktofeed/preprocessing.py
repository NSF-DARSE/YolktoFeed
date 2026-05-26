import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests


def significant_gene_filter(path: str,
                            omit_genes_without_symbol: bool = True,
                            min_val: float = 0.1,
                            min_count_per_day: int = 5,
                            omit_days: list = [],
                            use_cudf: bool = False,
                            verbose: int = 0
                            ):
    """
    Preprocess/filter the gene expression data based on minimum
    avergae, varaince across sample dataset
    
    Args:
        path (str): path to csv files folder
        omit_genes_without_symbol (bool):
        min_val (float):
        min_count_per_day (int):
        omit_days (list):
        use_cudf (bool): use gpu for data preprocessing

    Returns:
        pandas.DataFrame: filtered genes dataframe
        pandas.DataFrame: metadata df

    Raises:
        FileNotFoundError: if corresponding file can't be found in path
        
    Example:
        >>> gene_df = significant_gene_filter(path="data/")
    """
    
    df_base = pd.read_csv(path + 
                    "Liver_HEIDI_complete_ensembl_symbol.csv")

    if verbose: print("Number of genes: ", df_base.shape[0])

    # finding last row with valid gene symbol
    gene_list_last_index = (df_base["Gene Symbol"] == "-").idxmax()

    # printing number of nans
    if verbose: print("Number of nan's: ", 
                      df_base.iloc[:,2:].apply(pd.to_numeric, 
                        errors='coerce').isna().sum().sum())

    
    if omit_genes_without_symbol:
        # omitting genes without any names
        df_trunk = df_base.iloc[:gene_list_last_index]
        df_trunk = df_trunk.drop(columns=["Ensembl Gene ID"])
        df_trunk = df_trunk.set_index("Gene Symbol")
    else:
        df_trunk = df_base.copy()
        df_trunk.index = df_trunk['Gene Symbol'].\
            replace('-', pd.NA).fillna(df_trunk['Ensembl Gene ID'])

    
    #dropping some days
    drop_days=tuple(['_'+str(day) for day in omit_days])
    df_trunk = df_trunk.loc[:, ~df_trunk.columns.str.endswith(drop_days)]
    if verbose: print(f"df shape after dropping {omit_days} days:",
                                                 df_trunk.shape)
    
    
    meta_cols = ['Gene Symbol', 'Ensembl Gene ID']
    expr_cols = [c for c in df_trunk.columns if c not in meta_cols]

    sample_meta = pd.DataFrame({'Sample': expr_cols})
    sample_meta['Day'] = sample_meta['Sample'].str.extract(r'_(\d+)$').astype(int)

    #obtaining list of cols as pandas index object
    sample_cols = df_trunk.columns
    #obtaining days as an array
    days = sample_cols.str.extract(r'_(\d+)$')[0].values

    #getting a new df with all elemets satifying conditon
    X = df_trunk[sample_cols]
    expressed = X > min_val

    #getting genes that satisfy the condition grouped by day on
    #all days, replace with any for any day
    counts_by_day = expressed.T.groupby(days).sum().T
    keep = (counts_by_day >= min_count_per_day).all(axis=1)

    #keeping only those genes
    df_filtered = df_trunk.loc[keep]
    if verbose: print("Number of genes finale: ", df_filtered.shape[0])

    
    return df_filtered, sample_meta   

def pca_analysis(df_filtered):
    """
    runs PCA analysis on the filterd genes 
    Args:
        df_filtered (pandas.DataFrame):
        expr_cols (list[str]): 

    Returns:
        skleran.pca: 
        numpy.ndarray:

    Raises:
        FileNotFoundError: if corresponding file can't be found in path
        
    Example:
        >>> Xpca, pca_t = pca_analysis(df_filtered)
    
    """
    
    expr = df_filtered.copy()
    log_expr_vals = np.log2(expr.values + 1)
    log_expr = pd.DataFrame(log_expr_vals, index=expr.index, 
                                columns=expr.columns)
    
    X_t = log_expr.T.values.astype(float)
    X_t_scaled = StandardScaler().fit_transform(X_t)
    
    pca_t = PCA(n_components=10, random_state=42)
    Xpca  = pca_t.fit_transform(X_t_scaled)

    print(f"  PC1 explains: {pca_t.explained_variance_ratio_[0]*100:.1f}% of variance")
    print(f"  PC2 explains: {pca_t.explained_variance_ratio_[1]*100:.1f}% of variance")
    print(f"  PC3 explains: {pca_t.explained_variance_ratio_[2]*100:.1f}% of variance")
    print(f"  PC4 explains: {pca_t.explained_variance_ratio_[3]*100:.1f}% of variance")
    print(f"  Xpca shape  : {Xpca.shape}")

    return pca_t, Xpca

def plot_pca_components(pca_t, Xpca, sample_meta, DAYS: list):
    DAY_COLORS = {
        4:  '#2ecc71',   # green
        6:  '#27ae60',   # dark green
        8:  '#f1c40f',   # yellow
        10: '#e67e22',   # orange
        12: '#e74c3c',   # red
        14: '#8e44ad',   # PURPLE ← Day 14
        16: '#3498db',   # blue
        18: '#1abc9c',   # teal
        20: '#2c3e50',   # dark navy
        }

    fig, axes = plt.subplots(1, 1, figsize=(8, 8))
    fig.suptitle('PC1 vs PC2', fontsize=14, fontweight='bold')

    ax = axes
    for day in DAYS:
        idx = np.where(sample_meta['Day'].values == day)[0]
        ax.scatter(Xpca[idx, 0], Xpca[idx, 1],c=DAY_COLORS[day], 
                        label=f'Day {day}', s=90, 
                        edgecolors='white', linewidth=0.5, zorder=3)

        ax.set_xlabel(f'PC1 ({pca_t.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca_t.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_xlim((-300,300))
        ax.set_ylim((-100,100))
        ax.set_title(f'PCA for days in {DAYS}', fontsize=12)
        ax.legend(title='Day', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
        ax.grid(True, alpha=0.4)
    
    
    plt.tight_layout()
    plt.savefig(f"pca_components_in{DAYS}.png".replace(" ",""))
       


def differential_expression(df_filtered, sample_meta, padj_max=0.05, 
                                        days_a=[4], days_b=[8]):
    """
    Getting the differntial expression of the 
    """
    expr = df_filtered.copy()
    log_expr_vals = np.log2(expr.values + 1)
    log_expr = pd.DataFrame(log_expr_vals, index=expr.index, 
                                columns=expr.columns)
    
    day_a_samples = sample_meta[sample_meta['Day'].isin(days_a)]['Sample'].tolist()
    day_b_samples = sample_meta[sample_meta['Day'].isin(days_b)]['Sample'].tolist()

    #getting the mean across columns
    day_a_mean   = log_expr[day_a_samples].mean(axis=1)
    day_b_mean   = log_expr[day_b_samples].mean(axis=1)

    lfc_rest        = day_a_mean - day_b_mean
    prest_vals = []
    for g in log_expr.index:
        a = log_expr.loc[g, day_a_samples].values.astype(float)
        b = log_expr.loc[g, day_b_samples].values.astype(float)

        a = a[~np.isnan(a)]
        b = b[~np.isnan(b)]

        if len(a) < 2 or len(b) < 2:
            p = 1.0
        elif np.var(a) == 0 or np.var(b) == 0:
            p = 1.0 if np.mean(a) == np.mean(b) else 0.0
        else:
            res = stats.ttest_ind(a, b, equal_var=False)
            p = float(res.pvalue) if not np.isnan(res.pvalue) else 1.0

        prest_vals.append(p)

    
    de_df = pd.DataFrame({'Gene': log_expr.index, 
                        'LFC': lfc_rest.values, 'pval': prest_vals})
    #raw p values may produce false positives
    #benjamin-hochberg FDR false discover rate correction
    de_df['padj'] = multipletests(de_df['pval'], method='fdr_bh')[1]

    de_df['-log10p'] = -np.log10(np.clip(de_df['pval'], 1e-300, 1))
    de_df['sig']     = (de_df['padj'] < padj_max) & (de_df['LFC'].abs() > 1)

    de_df.attrs['days_a'] = days_a
    de_df.attrs['days_b'] = days_b

    up = de_df[de_df['sig'] & (de_df['LFC'] > 0)]
    dn = de_df[de_df['sig'] & (de_df['LFC'] < 0)]

    print(f"  Upregulated   : {len(up)} genes")
    print(f"  Downregulated: {len(dn)} genes")

    return de_df, up, dn

def volcano_plot(de_df, up, dn):
    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    ax = axes
    ns = de_df[~de_df['sig']]
    comp_days = f'{de_df.attrs['days_b']}vs{de_df.attrs['days_a']}'
    comp_days.replace(" ","")

    ax.scatter(ns['LFC'],  ns['-log10p'],  c='#bdc3c7', s=5,  alpha=0.4, label='Not significant')
    ax.scatter(up['LFC'],  up['-log10p'],  c='#e74c3c', s=15, alpha=0.8, label=f'Up at {de_df.attrs['days_b']} ({len(up)})')
    ax.scatter(dn['LFC'],  dn['-log10p'],  c='#3498db', s=15, alpha=0.8, label=f'Down at {de_df.attrs['days_b']} ({len(dn)})')

    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline( 1, color='#e74c3c', linestyle='--', linewidth=1, alpha=0.4)
    ax.axvline(-1, color='#3498db', linestyle='--', linewidth=1, alpha=0.4)

    for _, row in up.nlargest(4, '-log10p').iterrows():
        ax.annotate(row['Gene'], (row['LFC'], row['-log10p']),
             fontsize=7, color='#c0392b', xytext=(3, 2), textcoords='offset points')
    for _, row in dn.nlargest(4, '-log10p').iterrows():
        ax.annotate(row['Gene'], (row['LFC'], row['-log10p']),
             fontsize=7, color='#2980b9', xytext=(3, 2), textcoords='offset points')
    ax.set_xlabel(f'Log2 Fold Change ({de_df.attrs['days_b']} vs {de_df.attrs['days_a']})', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title('Volcano Plot', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.savefig(f"volcanoPlotFor{comp_days}.png")

path = "../data/"
gene_df, sample_meta = significant_gene_filter(path, 
                                omit_days=[6], min_val=0.1, 
                                min_count_per_day=5,verbose=1)

pca_t, Xpca = pca_analysis(gene_df)

plot_pca_components(pca_t, Xpca, sample_meta, DAYS=[4,8])

de_df, up, dn = differential_expression(gene_df, sample_meta, 
                                        padj_max=0.05, 
                                        days_a=[8], days_b=[10])
volcano_plot(de_df, up, dn)
print(up)


