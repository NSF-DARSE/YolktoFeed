def get_stats(df, col: str, limit=1e-6):
        """
        df is passed by reference, but the operation df[boolean]
        creates new dataframe
        
        Args:
            df (pandas.Dataframe)
            col (str): column to filter by
            limit (float): limiting value to filter by
        Returns:
            void
        """
        temp_df = df[df.loc[:,col] < limit]
        print("number of Genes with {} less than {}".format(col, limit), 
                temp_df.shape[0], 
                temp_df.loc[:,'sigma'].max(),
                temp_df.loc[:,'max'].max()) 

def significant_gene_filter(path: str,
                            min_avg = 0.1,
                            min_sigma = 0.1,
                            use_cudf = False,
                            verbose = 0
                            ):
    """
    Preprocess/filter the gene expression data based on minimum
    avergae, varaince across sample dataset
    
    Args:
        path (str): path to csv files folder
        min_avg (float): minimum averge to filter by
        min_sigma (float): minimum variance to filter by
        verbose (bool): print condition
        use_cudf (bool): use gpu for data preprocessing

    Returns:
        pandas.Dataframe: preprocessed data

    Raises:
        FileNotFoundError: if corresponding file can't be found in path
        
    Example:
        >>> gene_df = significant_gene_filter(path="data/")
    """
    if not use_cudf:
        import pandas as pd
    
    df_base = pd.read_csv(path + 
                    "Liver_HEIDI_complete_ensembl_symbol.csv")

    if verbose: print("Number of genes: ", df_base.shape[0])

    # finding last row with valid gene symbol
    gene_list_last_index = (df_base["Gene Symbol"] == "-").idxmax()

    # printing number of nans
    if verbose: print("Number of nan's: ", 
                      df_base.iloc[:,2:].apply(pd.to_numeric, 
                        errors='coerce').isna().sum().sum())

    # omitting genes without any names
    df_trunk = df_base.iloc[:gene_list_last_index]
    df_trunk = df_trunk.drop(columns=["Ensembl Gene ID"])
    df_trunk = df_trunk.set_index("Gene Symbol")

    
    # Adding a column for average gene expression across all sample 
    df_trunk.insert(0, "avg", df_trunk.iloc[:, 1:].sum(axis=1)\
                                    /(df_trunk.shape[1]-1))
    df_trunk.insert(1, "sigma", df_trunk.iloc[:,1:].var(axis=1))
    df_trunk.insert(2, "max", df_trunk.iloc[:,1:].max(axis=1))

    if verbose: get_stats(df_trunk, 'avg', 1e-6)
    if verbose: get_stats(df_trunk, 'avg', 1e-4)

    # filter based on avg, creates a new df    
    df = df_trunk[(df_trunk.loc[:,'avg'] > min_avg)]
    
    if verbose: get_stats(df, 'sigma', 1e-2)
    if verbose: get_stats(df, 'sigma', 1e-1)

    # filter based on sigma, creates a new df
    df = df[(df.loc[:,'sigma'] > min_sigma)]

    if verbose: print("Number of genes left: ", df.shape[0])

    # Some metrics of the final database
    if verbose: print("final, min of the avg values", df.loc[:,'avg'].min(),)
    if verbose: print("final, min sigma values", df.loc[:,'sigma'].min())
    
    return df    

path = "../data/"
gene_df = significant_gene_filter(path, verbose=1)
