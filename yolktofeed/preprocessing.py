import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt


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
        pandas.Dataframe: filtered genes dataframe
        list: list of expresssion cols

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

    
    return df_filtered, expr_cols, sample_meta   

def pca_analysis(df_filtered, expr_cols):
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
        >>> Xpca, pca_t = pca_analysis(df_filtered, expr_cols)
    
    """
    
    expr = df_filtered[expr_cols].copy()
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
       

path = "../data/"
gene_df, expr_cols, sample_meta = significant_gene_filter(path, 
                                    omit_days=[6], min_val=0.1, 
                                    min_count_per_day=5,verbose=1)

pca_t, Xpca = pca_analysis(gene_df, expr_cols)

plot_pca_components(pca_t, Xpca, sample_meta, DAYS=[4,8])
