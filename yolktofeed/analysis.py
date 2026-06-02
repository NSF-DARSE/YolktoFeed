import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def pca_analysis(df_filtered):
    """
    runs PCA analysis on the filterd genes 
    Args:
        df_filtered (pandas.DataFrame):

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


def differential_expression(df_filtered, sample_meta, padj_max=0.05, 
                                        days_a=[4], days_b=[8]):
    """
    Getting the differntial expression of the filtered genes of two
    pairs of days

    Args:
        df_filtered (pandas.DataFrame):
        sample_meta (pandas.DataFrame):
        padj_max (float):
        days_a (list[int]):
        days_b (list[int]):

    Returns:
        pandas.DataFrame: genes with de values
        pandas.DataFrame: up regulated genes
        pandas.DataFrame: dn regulated genes

    Raises:
        FileNotFoundError: if corresponding file can't be found in path
        
    Example:
        >>> de_df, up, dn = differntail_expression(df_filetered, 
        >>>                                  days_a=[4], days_b=[8])
    
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

        nan_a = np.isnan(a).sum()
        nan_b = np.isnan(b).sum()
        if nan_a > 0 or nan_b > 0:
            print(f"{g}: NaNs in a={nan_a}, b={nan_b}")
            
    
        #duplicate genes issues, this part would linearize array
        a = a[~np.isnan(a)]
        b = b[~np.isnan(b)]
        
        #sample size min of 2 required for test
        if len(a) < 2 or len(b) < 2:
            p = 1.0
            print("!Warningsample size less that 2:", g)
        #if varince is zero t-test is unstable
        elif np.var(a) == 0 or np.var(b) == 0:
            p = 1.0 if np.mean(a) == np.mean(b) else 0.0
            print("equal variance case: ", g)
        else:
            res = stats.ttest_ind(a, b, equal_var=False)
            p = float(res.pvalue) if not np.isnan(res.pvalue) else 1.0

        prest_vals.append(p)

    
    de_df = pd.DataFrame({'Gene': log_expr.index, 
                            'LFC': lfc_rest.values,
                            'absLFC':np.abs(lfc_rest.values), 
                            'pval': prest_vals})
    #raw p values may produce false positives
    #benjamin-hochberg FDR false discover rate correction
    de_df['padj'] = multipletests(de_df['pval'], method='fdr_bh')[1]

    de_df['-log10padj'] = -np.log10(np.clip(de_df['padj'], 1e-300, 1))
    de_df['sig']     = (de_df['padj'] < padj_max) & (de_df['LFC'].abs() > 1)

    de_df.attrs['days_a'] = days_a
    de_df.attrs['days_b'] = days_b

    up = de_df[de_df['sig'] & (de_df['LFC'] > 0)]
    dn = de_df[de_df['sig'] & (de_df['LFC'] < 0)]

    print(f"  Upregulated   : {len(up)} genes")
    print(f"  Downregulated: {len(dn)} genes")

    return de_df, up, dn

