# Compartment-associated expression in human diabetic nephropathy

Reproducible analysis code for the revised human diabetic nephropathy (DN)
spatial-transcriptomics manuscript. The primary analysis uses AOI-level
TMM–voom–dream mixed-effects models. A patient–compartment pseudobulk analysis
assesses sensitivity to aggregation within the same cohort.

- Data accession: [GSE346609](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE346609)
- Code repository: [Heee315/DN_Spatial](https://github.com/Heee315/DN_Spatial)

This release replaces the earlier DESeq2 workflow. It does not exclude the eight
previously flagged genes from the primary analysis.

## Files

```text
project_release/
├── README.md
├── VALIDATION.md
├── analysis.R
├── .gitignore
├── data/
│   ├── README.md
│   ├── raw_counts.csv                         # GEO input
│   ├── GSE346609_Slide_Assignment.xlsx        # GEO input
│   ├── input_checksums.csv
│   └── reference/
│       ├── README.md
│       ├── GO_displayed_terms.csv
│       ├── original_package_versions.csv
│       ├── original_sessionInfo.txt
│       └── PPI_*_mapping.csv / PPI_*_edges.csv
├── validation/                               # small release-check records
└── results/
    └── .gitkeep
```

The local release archive includes the two analysis inputs. In a Git checkout,
these inputs and generated results are ignored by `.gitignore`; obtain the
inputs from GEO and place them in `data/`. Commit the code, documentation and
small reference files. Upload the **contents** of `project_release/` to the
repository root when replacing the earlier workflow.

## Input data

`raw_counts.csv` must have one row per unique human gene symbol and one column per
AOI. The first column contains gene symbols; the remaining columns contain
unnormalized, nonnegative integer gene counts. This is the gene-level count
matrix used for analysis, not sequencing FASTQ data or normalized expression.

The slide-assignment workbook must contain these columns on its first sheet:

| Column | Use |
|---|---|
| `ROILabel` | AOI identifier; used to derive patient, disease and compartment |
| `ScanLabel` | Cross-check of disease assignment (`NC` or `DN`) |
| `SlideName` | Slide assignment and confounding audit |
| `AOISurfaceArea` | Recorded AOI area; retained as metadata |
| `AOINucleiCount` | Recorded nuclei count; retained as metadata |
| `ROICoordinateX`, `ROICoordinateY` | Recorded ROI coordinates |

Accepted names include `Control-1-G1`, `Normal-1-G1`,
`NC | Normal-1-G1`, `DN-2-G1` and `DN | DN-2-G1`.
`Normal` is standardized to `Control`. Counts and metadata are matched by AOI
identifier, never by their input order. Duplicates, unmatched AOIs, invalid counts
and unexpected cohort sizes stop the analysis.

The manuscript cohort comprises 46 AOIs from eight patients: 12 control G,
12 control T, 10 DN G and 12 DN T. DN-1 has no glomerular AOI. There are 15
patient–compartment profiles and seven patients with both compartments. No
missing glomerular profile is imputed.

See [data/README.md](data/README.md) for input provenance and QC scope.

## R packages

The original analyses used R 4.3.0. Recorded package versions are included in
`data/reference/original_package_versions.csv`. Important versions were
variancePartition 1.32.5, edgeR 3.42.4, limma 3.58.1, clusterProfiler 4.10.1 and
org.Hs.eg.db 3.18.0. This version list is a record of the original environment,
not a complete dependency lockfile. The original session report is also supplied
as `data/reference/original_sessionInfo.txt`, including annotation and model
dependency versions.

Install the required packages in R before running the script:

```r
install.packages(c("BiocManager", "readxl", "ggplot2", "ggrepel",
                   "patchwork", "igraph"))
BiocManager::install(c("edgeR", "limma", "variancePartition", "BiocParallel",
                       "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
                       "GEOquery"), ask = FALSE, update = FALSE)
# Only needed for an optional fresh STRING query:
install.packages("rbioapi")
```

Bioconductor installs a release compatible with the installed R version. A
current installation may differ from the original environment, especially in
annotation databases and model-fitting dependencies. For historical
reproduction, use R 4.3 and the recorded package versions from the corresponding
CRAN/Bioconductor archives. Do not interpret the installation block as an exact
historical environment restore.

The script records installed versions, compares available versions with the
original record (including model/optimizer dependencies), and audits selected
manuscript quantities. Differences are reported; results are never adjusted to
force agreement. Packages are not installed or updated by `analysis.R`.

The release validation reproduced the AOI interaction/attenuation gene lists and
all pseudobulk selected-gene lists. Some AOI fit failures and selected DEG counts
differed in the validation environment; see [VALIDATION.md](VALIDATION.md).

## Run

From a terminal in `project_release/`:

```bash
# First validate the two inputs without fitting any model:
Rscript analysis.R --check-only

# Run the main AOI and pseudobulk analyses, GO, cached STRING and figures:
Rscript analysis.R --workers=4
```

One worker is the default and is the lowest-memory option. Multiple workers use
BiocParallel SOCK workers. Fitting thousands of mixed-effects models can take
substantial time; optional sensitivity analyses add nine AOI refits.

After both supplementary inputs have been made public on GEO, they can be
downloaded and checked with:

```bash
Rscript analysis.R --download --check-only
Rscript analysis.R --workers=4
```

`--download` discovers supplementary filenames through GEOquery, including
GEO-prefixed names and a gzip-compressed count CSV. It requires public access to
the files; it cannot download a private deposition. To compare a newly downloaded
version with an existing local input, specify the new path explicitly.

Additional options:

```bash
Rscript analysis.R --help
Rscript analysis.R --counts=/path/raw_counts.csv --metadata=/path/GSE346609_Slide_Assignment.xlsx
Rscript analysis.R --workers=4 --output=results/main_run
Rscript analysis.R --workers=4 --sensitivity
Rscript analysis.R --workers=4 --skip-go --skip-ppi --skip-plots
Rscript analysis.R --workers=4 --string-live
```

Output directories must be new or empty. Without `--output`, each run creates a
timestamped directory in `results/`. Relative input/output overrides are
resolved from the working directory. The standard `data/` directory is resolved
relative to `analysis.R` when using `Rscript`.

## Analysis specification

### Primary AOI model and contrasts

All 46 AOIs are retained, including AOIs with QC flags. Genes with zero total
count are excluded as uninformative; no additional gene-level detection/LOQ
filter or predefined gene exclusion is applied. The supplied matrix contains
18,675 genes. TMM normalization precedes voomWithDreamWeights
(`prior.count=0.5`, `rescaleWeightsAfter=FALSE`, no additional between-array
normalization).

```r
~ 0 + Group + (1 | Patient_ID) + (1 | Patient_Compartment)
```

The four groups are control T, control G, DN T and DN G. Models use REML,
Satterthwaite degrees of freedom and empirical-Bayes moderation. The five direct
contrasts are:

| Output name | Contrast |
|---|---|
| `G_DN_vs_Control` | DN G − control G |
| `T_DN_vs_Control` | DN T − control T |
| `Control_GT` | Control G − control T |
| `DN_GT` | DN G − DN T |
| `Interaction` | (DN G − DN T) − (control G − control T) |

BH adjustment is applied separately to each contrast. Non-interaction DEGs
require FDR < 0.05 and absolute log2 fold change > 1. Interactions require
FDR < 0.05. Unavailable gene fits remain in full tables with missing statistics.
Confidence intervals are 95% pointwise moderated-t intervals calculated using
each gene/contrast's standard error and degrees of freedom.

Slide and disease are completely confounded in this cohort. A slide effect is
not fitted as if it were independently estimable. Slide-by-disease counts and
design rank are exported. A shared sequencing run does not resolve slide or
procurement confounding.

### Interaction and attenuation

Let `C` and `D` be control and DN G–T log2 fold changes, respectively.

| Class | Criteria |
|---|---|
| Glomerular attenuation | Significant interaction; control G–T FDR < 0.05; `C > 1`; `0 <= D < C` |
| Tubular attenuation | Significant interaction; control G–T FDR < 0.05; `C < -1`; `C < D <= 0` |
| Other interaction | Significant interaction outside either attenuation definition |
| Not significant | Interaction does not meet the significance criterion |

A remaining DN G–T contrast can be significant and still be attenuated.
Estimated sign reversals are annotated within the other-interaction class,
rather than made into an additional attenuation class.

### Pseudobulk, distances and markers

Raw counts are summed within each patient–compartment to generate 15 pseudobulk
profiles. TMM, voom and dream are refitted using:

```r
~ 0 + Group + (1 | Patient_ID)
```

AOI–pseudobulk comparisons use the same contrasts and gene identifiers.
Agreement evaluates aggregation within this cohort, not independent replication.

For Figure 2B, AOI logCPM values are averaged within patient–compartment. Each of
the control G–T DEGs is standardized across the 14 profiles from the seven
complete-pair patients. Each patient's G–T distance is the root mean square of
the gene-wise standardized differences. Group means are descriptive; no
confirmatory P value is generated for this same-cohort selected-gene summary.

For the 17 selected markers, both AOI and mean-profile heatmaps use each gene's
mean and sample SD across all 15 patient–compartment mean profiles. Heatmap colors
are capped at −3 and +3; exported z-scores are not capped. Marker expression
provides regional context and does not estimate cell fractions.

### GO and STRING

GO biological-process enrichment is performed separately for glomerular
attenuation, tubular attenuation, downregulated glomerular DN DEGs and
downregulated tubular DN DEGs. One-to-one SYMBOL–ENTREZ mappings are used;
ambiguous mappings are excluded from both foreground and background.
The background is the full testable gene set for the corresponding contrast,
with GO annotation defining the effective universe. Attenuation enrichment uses
the broad testable background, not only control compartment-associated genes.
This background can capture baseline compartment specialization.

GO testing uses gene-set sizes 10–500 and BH FDR < 0.05 within each analysis.
All tested terms and all significant terms are exported. The small reference
file `GO_displayed_terms.csv` controls only the terms displayed in Figures 3 and
S3; it does not select the genes or restrict the enrichment tests. Archived
display terms that are no longer significant are flagged and omitted.

By default, STRING networks are reconstructed from archived mapping/edge tables
to avoid changes in the live database. Input genes are selected from the current
AOI attenuation results. Human mappings, accepted identifiers, combined score
≥ 0.7 and zero additional nodes are enforced. The incorrect VEGFA→COL18A1 mapping
is rejected; G6PC1→G6PC is an accepted alias. Disconnected input proteins are
retained in node audit tables and omitted from the plotted networks.
If a changed analysis selects genes absent from the cache, the script stops
instead of silently dropping them. Use `--string-live` only for a deliberately
updated database query. Live results may differ.

STRING edges represent functional associations; the code does not infer causal
hubs, direct physical binding or disrupted interactions in the sampled kidneys.
See [data/reference/README.md](data/reference/README.md) for cache provenance.

## Outputs

| Location within a run | Contents |
|---|---|
| `AOI/` | Five complete contrast tables, selected DEGs and interaction classes |
| `Pseudobulk/` | Summed counts, fitted contrasts and interaction classes |
| `Aggregation/` | Matched effect estimates, threshold selections and Spearman correlations |
| `Distance/` | Patient RMS distances, summaries and gene standardization |
| `Markers/`, `Display/` | AOI and patient–compartment expression and z-scores |
| `GO/` | Mapping audit, foreground/background IDs and full/significant/displayed GO terms |
| `STRING/` | Mapping audit, node/edge tables and network summary |
| `Figures/` | PDF and 300-dpi PNG figures, PCA data and network coordinates |
| `Sensitivity/` | Optional complete-pair and leave-one-patient-out AOI fits |
| Run root | Settings, input checksums, versions, session information and manuscript comparison |

Figures are generated in manuscript order:

1. Figure 1: AOI PCA and control/DN G–T contrasts.
2. Figure 2: interaction classes, selected-gene RMS and representative attenuation contrasts.
3. Figure 3: compartment-specific DN–control DEGs and downregulated-DEG GO.
4. Figure 4: AOI marker heatmap.
5. Figure S2: AOI–pseudobulk effect comparisons.
6. Figure S3: attenuation-gene GO.
7. Figure S4: attenuation-gene STRING networks.
8. Figure S5: patient–compartment marker means and AOI expression displays.

Figure S1 contains microscopy images and is not reconstructible from these two
numeric input files. It remains a separately supplied manuscript figure.
Plots generated here reproduce the analysis and panel content; network layouts,
labels and final publication artwork need not be pixel-identical.

`manuscript_comparison.csv` checks selected expected results, including 291 AOI
interactions (145 G attenuation, 11 T attenuation and 135 other interactions),
335 pseudobulk interactions and the displayed RMS/correlation summaries.
`RUN_COMPLETE.txt` is written only after the requested stages finish. Review
warnings, missing fits, version comparisons and the manuscript audit before
using a new run for publication.

## Method documentation

- [variancePartition dream workflow](https://gabrielhoffman.github.io/variancePartition/articles/dream.html)
- [dream model reference](https://gabrielhoffman.github.io/variancePartition/reference/dream-method.html)
- [GEOquery](https://bioconductor.org/packages/GEOquery/)
- [STRING](https://string-db.org/)

This code does not reprocess sequencing reads, regenerate vendor probe-level QC,
perform cell-fraction deconvolution, or reconstruct tissue procurement history
from the count matrix.

See [VALIDATION.md](VALIDATION.md) for the checks performed on this release and
the distinction between fresh runs and checks using archived results.
