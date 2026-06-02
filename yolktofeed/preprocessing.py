import pandas as pd

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

"""

path = "../data/"

gene_df, sample_meta = significant_gene_filter(path, 
                                omit_days=[6], min_val=0.1, 
                                min_count_per_day=5,verbose=1)

pca_t, Xpca = pca_analysis(gene_df)

plot_pca_components(pca_t, Xpca, sample_meta, DAYS=[18,20])

de_df, up, dn = differential_expression(gene_df, sample_meta, 
                                        padj_max=0.05, 
                                        days_a=[18], days_b=[20])
volcano_plot(de_df, up, dn, sortBy='absLFC')

updnRegulatedGenes(de_df, up, dn, sortBy='absLFC')
#print(up)
"""
