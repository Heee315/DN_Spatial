# Release validation

This report distinguishes fresh computations from checks that reuse the
archived manuscript model estimates. Expected manuscript counts are audit
targets; they are never substituted for computed results.

## Input and normalization checks

- R syntax and the command-line help were checked.
- `--check-only` passed on the included 18,675-gene × 46-AOI matrix and workbook.
- AOI names/order and the supplied slide, area, nuclei and coordinate values
  agree with the archived analysis metadata after `Normal`→`Control` conversion.
- Reversing the count-column order and using prefixed `NC | Normal-*` / `DN | DN-*`
  identifiers returns exactly the same matched matrix.
- Duplicate AOIs and negative counts are rejected. Gzip-compressed CSV input
  returns the same count matrix as the uncompressed file.
- Raw-count pseudobulk aggregation produces 15 profiles and preserves per-gene
  total counts. No glomerular profile is imputed for DN-1.
- Independently recalculated TMM logCPM agrees with the archived voom expression;
  the largest absolute difference is approximately 5 × 10⁻¹⁴.

## Checks using archived model estimates

The release functions were executed on all archived AOI contrast tables and
the pseudobulk tables in the final supplementary workbook.

| Quantity | Recalculated value |
|---|---:|
| AOI significant interactions | 291 |
| AOI glomerular attenuation | 145 |
| AOI tubular attenuation | 11 |
| AOI other interactions, including estimated reversals | 135 |
| Pseudobulk significant interactions | 335 |
| Pseudobulk glomerular / tubular attenuation | 166 / 22 |
| Mean RMS, controls | 1.9605667743 |
| Mean RMS, DN | 0.9533502190 |
| AOI–pseudobulk Spearman rho, glomerular DN–control | 0.9876855988 |
| AOI–pseudobulk Spearman rho, tubular DN–control | 0.9943648234 |
| AOI–pseudobulk Spearman rho, interaction | 0.9792609093 |
| Cached glomerular STRING, connected proteins / edges | 75 / 113 |
| Cached tubular STRING, connected proteins / edges | 2 / 1 |

The GO figure specification was checked against the archived displayed-term
tables (15, 11, 15 and 15 terms). The new figure-generation functions produced
Figures 1–4 and S2–S5 from archived model/GO results, and the rendered figures
were visually inspected. These checks do not constitute a new GO enrichment
analysis or a new STRING database query.

## Fresh AOI and pseudobulk fits

The two model functions were also executed on the complete 18,675-gene input,
using eight workers. The same contrasts, thresholds and attenuation definitions
were applied; no reported value was forced into the results.

- The AOI interaction set (291 genes), glomerular-attenuation set (145) and
  tubular-attenuation set (11) match the archived gene lists exactly.
- All five pseudobulk selected-gene lists match the archived lists exactly:
  536 glomerular DN–control DEGs, 180 tubular DN–control DEGs, 488 control G–T
  DEGs, 97 DN G–T DEGs and 335 interactions. All 18,675 pseudobulk genes were testable.
- The cached STRING step was also rerun using the freshly selected AOI
  attenuation lists and reproduced 75/113 and 2/1 connected proteins/edges.

Some AOI model fits differed in availability in the validation environment:

| Quantity | Original run | Fresh validation run |
|---|---:|---:|
| AOI genes with available interaction tests | 18,657 | 18,654 |
| AOI glomerular DN–control DEGs | 480 | 481 |
| AOI tubular DN–control DEGs | 141 | 141 |
| AOI control G–T DEGs | 484 | 485 |
| AOI DN G–T DEGs | 58 | 58 |
| Mean control RMS | 1.9605667743 | 1.9613205395 |
| Mean DN RMS | 0.9533502190 | 0.9535131580 |

These selection changes involve genes whose fits were unavailable in one run
and available in the other, rather than simply crossing the FDR boundary.
The changed selections are listed in `validation/AOI_selection_changes.csv`.
Among genes fitted in both runs, the largest absolute AOI log2 fold-change
difference across the five contrasts was approximately 3.8 × 10⁻⁵. The RMS
summary uses each run's selected control G–T gene set, so it also changes when
that set changes. Failed fits remain missing; they are not replaced with
archived estimates or fitted with a simpler model.

The automated comparison and software record are supplied in `validation/`.
**Not every manuscript quantity was reproduced exactly in this environment.**
Use the recorded original R/package environment for exact historical
verification, and review model diagnostics and gene lists when versions differ.
The fresh validation outputs do not replace the manuscript's original results.

## Scope of the validation environment

The validation environment uses R 4.3.3 and variancePartition 1.32.5. Some
dependencies differ from the original R 4.3.0 environment; both the original
software record and per-run version comparisons are available.

Public GEO downloading, live STRING querying and the optional nine sensitivity
refits have not been exercised as part of release validation. Microscopy
Figure S1 requires its original image inputs. A fresh clusterProfiler enrichment
run has not been performed in this validation environment; the GO calls and
foreground/background definitions were checked against the original analysis
code and the clusterProfiler 4.10.1 source.
