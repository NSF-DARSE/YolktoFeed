import pandas as pd

path = ""
df_base = pd.read_csv(path+"Liver_HEIDI_complete_ensembl_symbol.csv")

print("Number of genes: ", df_base.shape[0])

# finding last row with valid gene symbol
gene_list_last_index = (df_base["Gene Symbol"] == "-").idxmax()

print("Number of nan's: ", df_base.iloc[:,2:].apply(pd.to_numeric, 
                                errors='coerce').isna().sum().sum())

#truncating till that row
df_trunk = df_base.iloc[:gene_list_last_index]
df_trunk = df_trunk.drop(columns=["Ensembl Gene ID"])
df_trunk = df_trunk.set_index("Gene Symbol")

# Adding a column for average gene expression across all days and all chikens for prelim analysis
df_trunk.insert(0, "avg", df_trunk.iloc[:, 1:].sum(axis=1)/(df_trunk.shape[1]-1))
df_trunk.insert(1, "sigma", df_trunk.iloc[:,1:].var(axis=1))
df_trunk.insert(2, "max", df_trunk.iloc[:,1:].max(axis=1))

df = df_trunk

Basic stats on average gene expression levels with their corresponding sigma
#column 0 is avg, 1 is sigma, 2 is max
ind = 0
print("number of Genes and max-variance, avg less than 1e-6 and its max sigma/val:\t", df.loc[df.iloc[:,ind] < 1e-6, df.columns[0:2]].shape[0], "\t",  df.loc[df.iloc[:,ind] < 1e-6, df.columns[2]].max(), "\t", df.loc[df.iloc[:,ind] < 1e-6, df.columns[3]].max())
print("number of Genes and max-variance, avg less than 1e-4 and its max sigma/val:\t", df.loc[df.iloc[:,ind] < 1e-4, df.columns[0:2]].shape[0], "\t",  df.loc[df.iloc[:,ind] < 1e-4, df.columns[2]].max(), "\t", df.loc[df.iloc[:,ind] < 1e-4, df.columns[3]].max())
print("number of Genes and max-variance, avg less than 1e-2 and its max sigma/val:\t", df.loc[df.iloc[:,ind] < 1e-2, df.columns[0:2]].shape[0], "\t",  df.loc[df.iloc[:,ind] < 1e-2, df.columns[2]].max(), "\t", df.loc[df.iloc[:,ind] < 1e-2, df.columns[3]].max())
print("number of Genes and max-variance, avg less than 1e-1 and its max sigma/val:\t", df.loc[df.iloc[:,ind] < 1e-1, df.columns[0:2]].shape[0], "\t",  df.loc[df.iloc[:,ind] < 1e-1, df.columns[2]].max(), "\t", df.loc[df.iloc[:,ind] < 1e-1, df.columns[3]].max())

df.loc[(df.loc[:,'avg'] < 1e-1) & (df.loc[:,'max'] > 5), df.columns[:3]]

df = df[(df.loc[:,'avg'] > 1e-1)]

print("Number of genes left: ", df.shape[0])

ind = 1
print("number of Genes with variance less than 1e-6:\t", df.loc[df.iloc[:,ind] < 1e-6, :].shape[0], "\t", df.loc[df.iloc[:,ind] < 1e-6, df.columns[3]].max())
print("number of Genes with variance less than 1e-4:\t", df.loc[df.iloc[:,ind] < 1e-4, :].shape[0], "\t", df.loc[df.iloc[:,ind] < 1e-4, df.columns[3]].max())
print("number of Genes with variance less than 1e-2:\t", df.loc[df.iloc[:,ind] < 1e-2, :].shape[0], "\t", df.loc[df.iloc[:,ind] < 1e-2, df.columns[3]].max())
print("number of Genes with variance less than 1e-1:\t", df.loc[df.iloc[:,ind] < 1e-1, :].shape[0], "\t", df.loc[df.iloc[:,ind] < 1e-1, df.columns[3]].max())

df = df[(df.loc[:,'sigma'] > 1e-1)]
print("Number of genes left: ", df.shape[0])

# Some metrics of the final database
print("Minimum of the avg values", df.loc[:,'avg'].min(),)
print("Minimum sigma values", df.loc[:,'sigma'].min())
#display(df.loc[df.loc[:,'sigma'] < 0.01, df.columns[:4]])
