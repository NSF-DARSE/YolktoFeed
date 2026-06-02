# YolktoFeed

Python package for preprocessing, PCA analysis, differential expression analysis, and volcano plot visualization.

## Create Virtual Environment

```bash
python -m venv .binfo
```

### Activate Environment

Linux/macOS:

```bash
source .binfo/bin/activate
```

Windows PowerShell:

```powershell
.binfo\Scripts\Activate.ps1
```

## Installation

```bash
pip install -e ".[dev]"
```

## Example Workflow

```python
from yolktofeed.preprocessing import significant_gene_filter
from yolktofeed.analysis import pca_analysis
from yolktofeed.analysis import differential_expression
from yolktofeed.plot2d import plot_pca_components, volcano_plot
from yolktofeed.plot2d import updnRegulatedGenes

path = </path/to/Liver_HEIDI_complete_ensembl_symbol.csv>

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
```

## Features

* Gene filtering and preprocessing
* PCA analysis and visualization
* Differential expression analysis
* Volcano plot generation
* Up/down-regulated gene summary tables

## Development

Run tests:

```bash
pytest
```

Build package:

```bash
python -m build
```
