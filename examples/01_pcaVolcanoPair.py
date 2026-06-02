from yolktofeed.preprocessing import significant_gene_filter
from yolktofeed.analysis import pca_analysis
from yolktofeed.analysis import differential_expression
from yolktofeed.plot2d import plot_pca_components, volcano_plot
from yolktofeed.plot2d import updnRegulatedGenes

path = "../data/"

gene_df, sample_meta = significant_gene_filter(
    path,
    omit_days=[6],
    min_val=0.1,
    min_count_per_day=5,
    verbose=1
)

pca_t, Xpca = pca_analysis(gene_df)

plot_pca_components(
    pca_t,
    Xpca,
    sample_meta,
    DAYS=[4, 8]
)

de_df, up, dn = differential_expression(
    gene_df,
    sample_meta,
    padj_max=0.05,
    days_a=[4],
    days_b=[8]
)
volcano_plot(de_df, up, dn, sortBy='absLFC')

updnRegulatedGenes(de_df, up, dn, sortBy='absLFC')

