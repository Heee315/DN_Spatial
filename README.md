# Spatial Transcriptomics Analysis

This repository contains the R code used for spatial transcriptomic analysis of glomerular and tubular compartments in diabetic nephropathy (DN).

## Data availability

The processed raw count matrix used in this analysis is available from the NCBI Gene Expression Omnibus (GEO) under accession number **GSE346609**.

Download `raw_counts.csv` from **GSE346609** and place it in the `data/` directory before running the analysis.

## Repository structure

```text
project_release/
├── README.md
├── analysis.R
├── data/
│   └── raw_counts.csv
└── results/
```

The `data/` and `results/` directories may contain `.gitkeep` files in the GitHub repository so that the empty directories are retained by Git.

## Input data

`raw_counts.csv` contains gene-level raw counts for **46 areas of interest (AOIs)** from glomerular and tubular compartments of control and DN kidney tissues.

The first column contains gene symbols, and the remaining columns contain AOI-level raw counts.

The dataset contains the following AOI groups:

- 12 control glomerular AOIs
- 12 control tubular AOIs
- 10 DN glomerular AOIs
- 12 DN tubular AOIs

Eight genes showing extreme within-group expression outliers suggestive of possible probe capture artifacts are retained in the GEO-deposited count matrix and excluded during downstream analysis as specified in `analysis.R`:

```text
FOXN2
GPX3
IGHG1
IGHG2
IGHG3
IGHG4
IGKC
IGLL5
```

## Requirements

The analysis was performed in R using the following packages:

```r
DESeq2
tidyverse
edgeR
ggVennDiagram
rbioapi
ggrepel
ggpubr
KEGGREST
AnnotationDbi
clusterProfiler
enrichplot
org.Hs.eg.db
eulerr
matrixStats
```

Install any missing packages before running the script.

## Running the analysis

Run the script from the root directory of the repository:

```bash
Rscript analysis.R
```

The script expects the input file at:

```text
data/raw_counts.csv
```

and writes generated output files to:

```text
results/
```

## Analyses included

The script performs the following analyses:

- quality-control-associated exclusion of predefined outlier genes
- principal component analysis
- relative log expression analysis
- renal microstructure marker expression analysis
- DESeq2 differential expression analysis
- volcano and MA plots
- analysis of loss of compartment-specific expression
- overlap analysis of glomerular and tubular DEGs
- Gene Ontology enrichment analysis
- KEGG enrichment analysis
- STRING protein–protein interaction analysis
- supplementary QC and expression plots

## Reproducibility

The script uses project-relative paths and does not depend on local server paths or institution-specific file locations.

Sample grouping is inferred from AOI column names rather than fixed column positions, allowing the analysis to remain reproducible even if the order of columns in `raw_counts.csv` changes.

The analysis script also saves session information to:

```text
results/sessionInfo.txt
```

for reproducibility.

## Data source

NCBI Gene Expression Omnibus (GEO): **GSE346609**
