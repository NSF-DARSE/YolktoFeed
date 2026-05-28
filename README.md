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
from yolktofeed.preprocessing import pca_analysis, plot_pca_components
from yolktofeed.preprocessing import differential_expression
from yolktofeed.plot2d import volcano_plot
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
    DAYS=[18, 20]
)

de_df, up, dn = differential_expression(
    gene_df,
    sample_meta,
    padj_max=0.05,
    days_a=[18],
    days_b=[20]
)

volcano_plot(de_df, up, dn)

updnRegulatedGenes(de_df, up, dn)
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
