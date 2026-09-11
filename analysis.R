# ============================================================
# Spatial transcriptomics analysis
# Public reproducibility script
#
# Expected repository structure:
#
# project_release/
# ├── analysis.R
# ├── data/
# │   └── raw_counts.csv
# └── results/
#
# The script uses only project-relative paths. No local server,
# institution, user, or original report paths are required.
#
# Input:
#   data/raw_counts.csv
#   - first column: gene symbols
#   - remaining 46 columns: AOI raw counts
#
# Analysis:
#   - QC-associated exclusion of 8 predefined outlier genes
#   - DESeq2 differential expression
#   - PCA and expression plots
#   - compartment-specific expression-loss analysis
#   - GO/KEGG enrichment
#   - STRING PPI networks
#   - supplementary QC plots
#
# NOTE:
# This script preserves the AOI-level statistical design used for
# the manuscript analysis so that the reported results are reproducible.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(DESeq2)
library(tidyverse)
library(edgeR)
library(ggVennDiagram)
library(rbioapi)
library(ggrepel)
library(ggpubr)
library(KEGGREST)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(eulerr)


# ============================================================
# 0-1. Unified group colors
#
# Use the same group colors throughout all plots in which
# color/fill represents Control_G, Control_T, DN_G, or DN_T.
# ============================================================

group_colors <- c(
  "Control_G" = "green4",
  "Control_T" = "slateblue2",
  "DN_G" = "brown2",
  "DN_T" = "goldenrod"
)


# ============================================================
# 1. Project directories and input file
# ============================================================

project_dir <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

data_dir <- file.path(
  project_dir,
  "data"
)

output_dir <- file.path(
  project_dir,
  "results"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

raw_count_file <- file.path(
  data_dir,
  "raw_counts.csv"
)

if (!file.exists(raw_count_file)) {
  stop(
    "Input file was not found: ",
    raw_count_file,
    "\nRun this script from the repository root containing data/raw_counts.csv."
  )
}

# Keep legacy output calls simple by writing all generated files
# into the results directory.
setwd(output_dir)


# ============================================================
# 2. Read GEO processed raw-count matrix
# ============================================================

raw_counts <- read.csv(
  raw_count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(raw_counts) < 2) {
  stop("raw_counts.csv must contain one gene column and AOI count columns.")
}

colnames(raw_counts)[1] <- "Gene"

if (anyNA(raw_counts$Gene) || any(raw_counts$Gene == "")) {
  stop("Missing gene symbols were found in raw_counts.csv.")
}

if (anyDuplicated(raw_counts$Gene) > 0) {
  stop("Duplicated gene symbols were found in raw_counts.csv.")
}

gene_names <- raw_counts$Gene

id <- raw_counts[
  ,
  -1,
  drop = FALSE
]

id[] <- lapply(
  id,
  function(x) {
    suppressWarnings(
      as.numeric(x)
    )
  }
)

id <- as.matrix(id)

if (anyNA(id)) {
  stop("NA values were generated while converting count columns to numeric.")
}

if (any(id < 0)) {
  stop("Negative counts were found in raw_counts.csv.")
}

if (any(id %% 1 != 0)) {
  stop("Non-integer values were found in raw_counts.csv.")
}

storage.mode(id) <- "integer"

rownames(id) <- gene_names

rm(
  raw_counts,
  gene_names
)

cat(
  "Raw count matrix:",
  nrow(id),
  "genes x",
  ncol(id),
  "AOIs\n"
)

if (ncol(id) != 46) {
  warning(
    "The manuscript analysis contains 46 AOIs, but the input matrix contains ",
    ncol(id),
    " AOIs."
  )
}


# ============================================================
# 3. Build AOI metadata from public sample names
#
# Expected examples:
#   DN | DN-2-G1
#   NC | Normal-1-G1
#   DN | DN-1-T1
#   NC | Normal-1-T1
#
# Sample order is NOT assumed.
# ============================================================

metaData <- data.frame(
  id = colnames(id),
  stringsAsFactors = FALSE
)

metaData$ShortLabel <- stringr::str_trim(
  stringr::str_remove(
    metaData$id,
    "^(DN|NC)\\s*\\|\\s*"
  )
)

metaData$Condition <- dplyr::case_when(
  stringr::str_detect(
    metaData$id,
    "^DN\\s*\\|"
  ) ~ "DN",
  stringr::str_detect(
    metaData$id,
    "^NC\\s*\\|"
  ) ~ "Control",
  TRUE ~ NA_character_
)

metaData$Compartment <- dplyr::case_when(
  stringr::str_detect(
    metaData$ShortLabel,
    "-G[0-9]+$"
  ) ~ "G",
  stringr::str_detect(
    metaData$ShortLabel,
    "-T[0-9]+$"
  ) ~ "T",
  TRUE ~ NA_character_
)

metaData$sample <- dplyr::case_when(
  metaData$Condition == "Control" &
    metaData$Compartment == "G" ~ "Control_G",

  metaData$Condition == "Control" &
    metaData$Compartment == "T" ~ "Control_T",

  metaData$Condition == "DN" &
    metaData$Compartment == "G" ~ "DN_G",

  metaData$Condition == "DN" &
    metaData$Compartment == "T" ~ "DN_T",

  TRUE ~ NA_character_
)

if (
  anyNA(metaData$Condition) ||
    anyNA(metaData$Compartment) ||
    anyNA(metaData$sample)
) {
  bad_ids <- metaData$id[
    is.na(metaData$Condition) |
      is.na(metaData$Compartment) |
      is.na(metaData$sample)
  ]

  stop(
    "Some AOI names could not be parsed:\n",
    paste(
      bad_ids,
      collapse = "\n"
    )
  )
}

metaData$sample <- factor(
  metaData$sample,
  levels = c(
    "Control_G",
    "Control_T",
    "DN_G",
    "DN_T"
  )
)

rownames(metaData) <- metaData$id

stopifnot(
  identical(
    colnames(id),
    rownames(metaData)
  )
)

print(
  table(metaData$sample)
)

expected_group_counts <- c(
  Control_G = 12,
  Control_T = 12,
  DN_G = 10,
  DN_T = 12
)

observed_group_counts <- table(metaData$sample)

if (
  !identical(
    as.integer(observed_group_counts[names(expected_group_counts)]),
    as.integer(expected_group_counts)
  )
) {
  warning(
    "AOI group counts differ from the manuscript dataset."
  )
}


# ============================================================
# 4. Preserve the GEO matrix before analysis-specific exclusion
# ============================================================

id_unfiltered <- id

out <- c(
  "FOXN2",
  "GPX3",
  "IGHG1",
  "IGHG2",
  "IGHG3",
  "IGHG4",
  "IGKC",
  "IGLL5"
)

missing_outlier_genes <- setdiff(
  out,
  rownames(id)
)

if (length(missing_outlier_genes) > 0) {
  warning(
    "Expected QC outlier genes missing from input: ",
    paste(
      missing_outlier_genes,
      collapse = ", "
    )
  )
}

id <- id[
  !(rownames(id) %in% out),
  ,
  drop = FALSE
]

cat(
  "Genes after QC-associated outlier exclusion:",
  nrow(id),
  "\n"
)

if (nrow(id) != 18667) {
  warning(
    "The manuscript reports 18,667 genes after exclusion, but the current matrix contains ",
    nrow(id),
    " genes."
  )
}


# ============================================================
# 5. Compartment-specific matrices
#
# Selection is metadata-based so results do not depend on
# column order in raw_counts.csv.
# ============================================================

id_G <- id[
  ,
  metaData$Compartment == "G",
  drop = FALSE
]

id_T <- id[
  ,
  metaData$Compartment == "T",
  drop = FALSE
]

G_meta <- metaData[
  colnames(id_G),
  ,
  drop = FALSE
]

T_meta <- metaData[
  colnames(id_T),
  ,
  drop = FALSE
]

stopifnot(
  ncol(id_G) == 22,
  ncol(id_T) == 24,
  identical(
    colnames(id_G),
    rownames(G_meta)
  ),
  identical(
    colnames(id_T),
    rownames(T_meta)
  )
)


# ============================================================
# 6. DESeq2 run
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = id,
  colData = metaData,
  design = ~ sample
)

dds <- DESeq(
  dds,
  minReplicatesForReplace = Inf
)


# ============================================================
# 7. PCA
#    Keep only PCA version 2, as previously selected
# ============================================================

vsdata <- vst(
  dds,
  blind = FALSE
)

pca_data <- plotPCA(
  vsdata,
  intgroup = "sample",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_data, "percentVar")
)

pca_data$ShortLabel <- metaData[
  pca_data$name,
  "ShortLabel"
]

pca_data$sample <- factor(
  pca_data$sample,
  levels = c(
    "Control_G",
    "Control_T",
    "DN_G",
    "DN_T"
  )
)

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    shape = sample,
    color = sample,
    label = ShortLabel
  )
) +
  geom_point(size = 3) +
  geom_text(
    vjust = -1,
    size = 3,
    show.legend = FALSE
  ) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "%"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "%"
    )
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5)
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  scale_shape_manual(
    values = c(
      "Control_G" = 15,
      "Control_T" = 15,
      "DN_G" = 17,
      "DN_T" = 17
    )
  )

ggsave(
  "02_PCA_plot_shape.png",
  plot = p_pca,
  width = 7,
  height = 5,
  bg = "white"
)


# ============================================================
# 8. RLE
#    Kept because user requested remaining analyses to remain
# ============================================================

vsdata_rle <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(vsdata_rle)

# Explicit namespace fixes rowMedians-not-found errors.
rle <- vst_matrix - matrixStats::rowMedians(vst_matrix)

group <- metaData$sample

png(
  "03_RLE_boxplot.png",
  width = 1100,
  height = 450
)

par(
  mar = c(10, 5, 4, 2)
)

boxplot(
  rle,
  outline = FALSE,
  ylab = "Relative Log Expression",
  las = 2,
  col = group_colors[as.character(group)],
  cex.axis = 1.3,
  cex.lab = 1.7,
  ylim = c(-2.5, 2.5),
  yaxt = "n"
)

axis(
  2,
  cex.axis = 1.5
)

legend(
  "topright",
  legend = names(group_colors),
  fill = group_colors,
  bty = "n",
  title = "",
  cex = 1.4,
  horiz = TRUE,
  inset = c(0, -0.5),
  xpd = TRUE,
  text.width = 8.3,
  x.intersp = 0.1,
  y.intersp = 0.8
)

dev.off()


# ============================================================
# 9. CPM helper functions
#
# DGEobj.utils::convertCounts() was replaced because it produced
# a cpm() namespace/method error in the current environment.
# edgeR::cpm() is called explicitly instead.
# ============================================================

calculate_cpm_log <- function(data) {

  data <- as.matrix(data)
  storage.mode(data) <- "numeric"

  if (
    nrow(data) == 0 ||
      ncol(data) == 0
  ) {
    stop(
      "calculate_cpm_log(): input matrix has zero rows or columns."
    )
  }

  edgeR::cpm(
    data,
    log = TRUE,
    prior.count = 0.25
  )
}


calculate_cpm <- function(data) {

  data <- as.matrix(data)
  storage.mode(data) <- "numeric"

  if (
    nrow(data) == 0 ||
      ncol(data) == 0
  ) {
    stop(
      "calculate_cpm(): input matrix has zero rows or columns."
    )
  }

  edgeR::cpm(
    data,
    log = FALSE
  )
}


# ============================================================
# 10. Log2 CPM
# ============================================================

dn_index <- which(
  G_meta$sample == "DN_G"
)

nc_index <- which(
  G_meta$sample == "Control_G"
)

dt_index <- which(
  T_meta$sample == "DN_T"
)

nt_index <- which(
  T_meta$sample == "Control_T"
)

stopifnot(
  length(dn_index) == 10,
  length(nc_index) == 12,
  length(dt_index) == 12,
  length(nt_index) == 12
)

dn_cols_log <- calculate_cpm_log(
  id_G[, dn_index, drop = FALSE]
)

nc_cols_log <- calculate_cpm_log(
  id_G[, nc_index, drop = FALSE]
)

dt_cols_log <- calculate_cpm_log(
  id_T[, dt_index, drop = FALSE]
)

nt_cols_log <- calculate_cpm_log(
  id_T[, nt_index, drop = FALSE]
)


# Convert each group to a long table.
# The original code split sample names with names_sep = " \\| ".
# Because sample names can contain multiple separators, use a safer
# one-column sample-name approach and assign the group explicitly.

make_expression_long <- function(
    expression_matrix,
    group_name
) {

  as.data.frame(
    expression_matrix
  ) %>%
    rownames_to_column(
      "Gene"
    ) %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Expression"
    ) %>%
    mutate(
      Group = group_name
    )
}


dn_long <- make_expression_long(
  dn_cols_log,
  "DN_G"
)

nc_long <- make_expression_long(
  nc_cols_log,
  "Control_G"
)

dt_long <- make_expression_long(
  dt_cols_log,
  "DN_T"
)

nt_long <- make_expression_long(
  nt_cols_log,
  "Control_T"
)

expression_long <- bind_rows(
  nc_long,
  nt_long,
  dn_long,
  dt_long
)

expression_long$Group <- factor(
  expression_long$Group,
  levels = c(
    "Control_G",
    "Control_T",
    "DN_G",
    "DN_T"
  )
)


# ============================================================
# 11. Renal marker expression
#
# FIX:
# Original code selected IGHA1/JCHAIN/... but then set factor
# levels to PODXL/NPHS1/... causing the marker Gene column to
# become NA. The selection is now matched to the intended renal
# marker panel used in Figure S1.
# ============================================================

renal_marker_genes <- c(
  "PODXL",
  "NPHS1",
  "NPHS2",
  "LRP2",
  "SLC34A1",
  "SLC9A3",
  "CALB1",
  "AQP2"
)

marker <- expression_long %>%
  filter(
    Gene %in% renal_marker_genes
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = renal_marker_genes
    )
  )


# ============================================================
# 12. Marker Boxplot
#     KEPT as requested
# ============================================================

p_marker_boxplot <- ggplot(
  marker,
  aes(
    x = Gene,
    y = Expression,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.6,
    position = position_dodge(
      width = 0.75
    )
  ) +
  geom_jitter(
    aes(
      color = Group
    ),
    size = 1.2,
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    )
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    x = NULL,
    y = "Log2 Expression (CPM)"
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 13,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      color = "black",
      size = 13
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12
    )
  )

p_marker_boxplot

ggsave(
  file.path(
    "04_Marker_boxplot.png"
  ),
  plot = p_marker_boxplot,
  width = 16,
  height = 6,
  bg = "white"
)

# ============================================================
# 13. Marker Violin Plot
#     KEPT and corrected
# ============================================================

p_marker_violin <- ggplot(
  marker,
  aes(
    x = Gene,
    y = Expression,
    fill = Group
  )
) +
  geom_violin(
    width = 0.70,
    trim = TRUE,
    scale = "width",
    alpha = 0.8,
    color = "gray30",
    position = position_dodge(
      width = 0.80
    )
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  labs(
    x = NULL,
    y = "Log2 Expression (CPM)",
    fill = NULL
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 13,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      color = "black",
      size = 13
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 12
    )
  )

p_marker_violin

ggsave(
  file.path(
    "05_Marker_violin.png"
  ),
  plot = p_marker_violin,
  width = 16,
  height = 6,
  bg = "white"
)

p_marker_violin_box <- ggplot(
  marker,
  aes(
    x = Gene,
    y = Expression,
    fill = Group
  )
) +
  
  # 1. Violin plot
  geom_violin(
    aes(
      group = interaction(
        Gene,
        Group
      )
    ),
    width = 0.8,
    alpha = 0.45,
    trim = TRUE,
    scale = "width",
    color = NA,
    position = position_dodge(
      width = 0.85
    )
  ) +
  
  # 2. Boxplot
  geom_boxplot(
    aes(
      group = interaction(
        Gene,
        Group
      )
    ),
    width = 0.15,
    alpha = 0.75,
    outlier.shape = NA,
    position = position_dodge(
      width = 0.85
    )
  ) +
  
  # 3. Individual sample points
  geom_point(
    aes(
      color = "black",
      group = interaction(
        Gene,
        Group
      )
    ),
    size = 1.5,
    alpha = 0,
    position = position_jitterdodge(
      jitter.width = 0.05,
      dodge.width = 0.85
    )
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_color_manual(
    values = group_colors
  ) +
  
  labs(
    x = NULL,
    y = "Log2 Expression (CPM)",
    fill = NULL,
    color = NULL
  ) +
  
  theme_classic(
    base_size = 15
  ) +
  
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 13,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      color = "black",
      size = 13
    ),
    axis.title.y = element_text(
      size = 15
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 12
    )
  )
ggsave(
  file.path(
    "03_Marker_violin_boxplot.png"
  ),
  plot = p_marker_violin_box,
  width = 20,
  height = 6,
  bg = "white"
)
# ============================================================
# 14. Pattern Boxplot
#     KEPT as requested
# ============================================================

median_df <- marker %>%
  group_by(
    Gene,
    Group
  ) %>%
  summarise(
    median_val = median(Expression),
    .groups = "drop"
  )

p_pattern_box <- ggplot(
  marker,
  aes(
    x = Group,
    y = Expression,
    fill = Gene
  )
) +
  geom_boxplot(
    position = "dodge",
    outlier.size = 1.5,
    outlier.alpha = 0.7,
    width = 0.7
  ) +
  geom_line(
    data = median_df,
    aes(
      x = Group,
      y = median_val,
      group = Gene,
      color = Gene
    ),
    position = position_dodge(width = 0.7),
    linewidth = 1.2
  ) +
  geom_point(
    data = median_df,
    aes(
      x = Group,
      y = median_val,
      color = Gene
    ),
    position = position_dodge(width = 0.7),
    size = 2
  ) +
  labs(
    x = " ",
    y = "Log2 Expression (CPM)",
    title = ""
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    axis.title = element_text(size = 17),
    axis.text = element_text(
      size = 15,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  )

ggsave(
  "06_Pattern_boxplot.png",
  plot = p_pattern_box,
  width = 8,
  height = 5,
  bg = "white"
)


# ============================================================
# 15. DEG analyses
#
# DR-positive/negative analysis has been completely removed.
# ============================================================

res_NC_DN_G <- results(
  dds,
  contrast = c(
    "sample",
    "DN_G",
    "Control_G"
  ),
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)

res_NC_DN_T <- results(
  dds,
  contrast = c(
    "sample",
    "DN_T",
    "Control_T"
  ),
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)

res_NC_NC <- results(
  dds,
  contrast = c(
    "sample",
    "Control_G",
    "Control_T"
  ),
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)

res_DN_DN <- results(
  dds,
  contrast = c(
    "sample",
    "DN_G",
    "DN_T"
  ),
  cooksCutoff = FALSE,
  independentFiltering = FALSE
)


# ============================================================
# 16. Helper: DESeq2 result table
# ============================================================

result_to_df <- function(
    res,
    comparison
) {

  out_df <- as.data.frame(
    res
  ) %>%
    rownames_to_column(
      "gene"
    ) %>%
    mutate(
      comparison = comparison,
      diff = case_when(
        !is.na(padj) &
          padj < 0.01 &
          log2FoldChange > 1 ~ "UP",

        !is.na(padj) &
          padj < 0.01 &
          log2FoldChange < -1 ~ "Down",

        TRUE ~ "None"
      )
    )

  out_df
}


DEG_G <- result_to_df(
  res_NC_DN_G,
  "DN_G_vs_Control_G"
)

DEG_T <- result_to_df(
  res_NC_DN_T,
  "DN_T_vs_Control_T"
)

DEG_Control_GT <- result_to_df(
  res_NC_NC,
  "Control_G_vs_Control_T"
)

DEG_DN_GT <- result_to_df(
  res_DN_DN,
  "DN_G_vs_DN_T"
)


# ============================================================
# 17. Save complete DESeq2 result tables
# ============================================================

write.csv(
  DEG_G,
  "DEG_DN_G_vs_Control_G_all_genes.csv",
  row.names = FALSE
)

write.csv(
  DEG_T,
  "DEG_DN_T_vs_Control_T_all_genes.csv",
  row.names = FALSE
)

write.csv(
  DEG_Control_GT,
  "DEG_Control_G_vs_Control_T_all_genes.csv",
  row.names = FALSE
)

write.csv(
  DEG_DN_GT,
  "DEG_DN_G_vs_DN_T_all_genes.csv",
  row.names = FALSE
)

all_four_contrast_results <- bind_rows(
  DEG_G,
  DEG_T,
  DEG_Control_GT,
  DEG_DN_GT
)

write.csv(
  all_four_contrast_results,
  "All_gene_DESeq2_results_4_contrasts.csv",
  row.names = FALSE
)


# ============================================================
# 18. Table S2
#
# Requested change:
# Do NOT restrict Table S2 to top 10 genes.
#
# File 1 contains EVERY tested gene from:
#   DN_G vs Control_G
#   DN_T vs Control_T
#
# File 2 contains only significant DEGs from those two contrasts,
# but includes ALL significant DEGs rather than top 10.
# ============================================================

Table_S2_all_genes <- bind_rows(
  DEG_G,
  DEG_T
) %>%
  select(
    comparison,
    gene,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj,
    diff
  ) %>%
  arrange(
    comparison,
    diff,
    padj
  )

write.csv(
  Table_S2_all_genes,
  "Table_S2_all_genes_DN_vs_Control.csv",
  row.names = FALSE
)

Table_S2_all_significant_DEGs <- Table_S2_all_genes %>%
  filter(
    diff != "None"
  )

write.csv(
  Table_S2_all_significant_DEGs,
  "Table_S2_all_significant_DEGs.csv",
  row.names = FALSE
)


# ============================================================
# 19. Volcano plot function
# ============================================================

make_volcano_plot <- function(
    deg_df,
    title,
    highlight_genes = character(0),
    up_label = "UP",
    down_label = "Down",
    up_color,
    down_color
) {

  plot_df <- deg_df %>%
    mutate(
      Volcano_group = case_when(
        diff == "UP" ~ up_label,
        diff == "Down" ~ down_label,
        TRUE ~ "None"
      )
    )

  label_df <- plot_df %>%
    filter(
      gene %in% highlight_genes
    )

  label_none <- label_df %>%
    filter(
      Volcano_group == "None"
    )

  label_sig <- label_df %>%
    filter(
      Volcano_group != "None"
    )

  ggplot(
    plot_df,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      color = Volcano_group
    )
  ) +
    geom_point(
      alpha = 0.75,
      size = 1.5,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = setNames(
        c(
          "lightgray",
          up_color,
          down_color
        ),
        c(
          "None",
          up_label,
          down_label
        )
      )
    ) +
    geom_text_repel(
      data = label_sig,
      aes(
        label = gene
      ),
      size = 4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    # Non-significant gene labels
    geom_text_repel(
      data = label_none,
      aes(
        label = gene
      ),
      color = "grey55",
      size = 4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    labs(
      title = title,
      x = "log2 fold change",
      y = "-log10(adj p-value)",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text.x = element_text(
        color = "black",
        size = 15
      ),
      axis.text.y = element_text(
        color = "black",
        size = 15
      ),
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      legend.text = element_text(size = 14)
    )
}



# ============================================================
# 20. MA plot function
# ============================================================

make_ma_plot <- function(
    deg_df,
    title,
    up_label = "UP",
    down_label = "Down",
    highlight_genes = character(0),
    up_color,
    down_color
) {

  plot_df <- deg_df %>%
    mutate(
      MA_group = case_when(
        diff == "UP" ~ up_label,
        diff == "Down" ~ down_label,
        TRUE ~ "None"
      )
    )

  label_df <- plot_df %>%
    filter(
      gene %in% highlight_genes
    )

  ggplot(
    plot_df,
    aes(
      x = log2(baseMean + 1),
      y = log2FoldChange,
      color = MA_group
    )
  ) +
    geom_point(
      alpha = 0.7,
      size = 1.8,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = setNames(
        c(
          "lightgray",
          up_color,
          down_color
        ),
        c(
          "None",
          up_label,
          down_label
        )
      )
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 4,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = "grey75",
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    labs(
      title = title,
      x = "Average expression (log2)",
      y = "log2 Fold Change",
      color = NULL
    ) +
    theme_classic(
      base_size = 15
    ) +
    theme(
      axis.title = element_text(size = 18),
      axis.text.x = element_text(
        color = "black",
        size = 15
      ),
      axis.text.y = element_text(
        color = "black",
        size = 15
      ),
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      legend.text = element_text(size = 14)
    )
}


# ============================================================
# 21. Volcano plots for all four manuscript contrasts
# ============================================================

candidate_genes <- c(
  "FN1",
  "IGHA1",
  "PIGR",
  "UMOD"
)

p_vol_G <- make_volcano_plot(
  DEG_G,
  "DN_G vs Control_G",
  up_label = "UP in DN_G",
  down_label = "UP in Control_G",
  highlight_genes = candidate_genes,
  up_color = group_colors[["DN_G"]],
  down_color = group_colors[["Control_G"]]
)

p_vol_T <- make_volcano_plot(
  DEG_T,
  "DN_T vs Control_T",
  up_label = "UP in DN_T",
  down_label = "UP in Control_T",
  highlight_genes = candidate_genes,
  up_color = group_colors[["DN_T"]],
  down_color = group_colors[["Control_T"]]
)

p_vol_Control_GT <- make_volcano_plot(
  DEG_Control_GT,
  "Control_G vs Control_T",
  up_label = "UP in Control_G",
  down_label = "UP in Control_T",
  highlight_genes = candidate_genes,
  up_color = group_colors[["Control_G"]],
  down_color = group_colors[["Control_T"]]
)

p_vol_DN_GT <- make_volcano_plot(
  DEG_DN_GT,
  "DN_G vs DN_T",
  up_label = "UP in DN_G",
  down_label = "UP in DN_T",
  highlight_genes = candidate_genes,
  up_color = group_colors[["DN_G"]],
  down_color = group_colors[["DN_T"]]
)

ggsave(
  "07_Volcano_DN_G_vs_Control_G.png",
  plot = p_vol_G,
  width = 7,
  height = 6,
  bg = "white"
)

ggsave(
  "08_Volcano_DN_T_vs_Control_T.png",
  plot = p_vol_T,
  width = 7,
  height = 6,
  bg = "white"
)

ggsave(
  "09_Volcano_Control_G_vs_Control_T.png",
  plot = p_vol_Control_GT,
  width = 7,
  height = 6,
  bg = "white"
)

ggsave(
  "10_Volcano_DN_G_vs_DN_T.png",
  plot = p_vol_DN_GT,
  width = 7,
  height = 6,
  bg = "white"
)


# ============================================================
# 22. MA plots for all four manuscript contrasts
# ============================================================

p_ma_G <- make_ma_plot(
  DEG_G,
  "DN_G vs Control_G",
  up_label = "UP",
  down_label = "Down",
  highlight_genes = candidate_genes,
  up_color = group_colors[["DN_G"]],
  down_color = group_colors[["Control_G"]]
)

p_ma_T <- make_ma_plot(
  DEG_T,
  "DN_T vs Control_T",
  up_label = "UP",
  down_label = "Down",
  highlight_genes = candidate_genes,
  up_color = group_colors[["DN_T"]],
  down_color = group_colors[["Control_T"]]
)

p_ma_Control_GT <- make_ma_plot(
  DEG_Control_GT,
  "Control_G vs Control_T",
  up_label = "UP in Control_G",
  down_label = "UP in Control_T",
  highlight_genes = candidate_genes,
  up_color = group_colors[["Control_G"]],
  down_color = group_colors[["Control_T"]]
)

p_ma_DN_GT <- make_ma_plot(
  DEG_DN_GT,
  "DN_G vs DN_T",
  up_label = "UP in DN_G",
  down_label = "UP in DN_T",
  highlight_genes = candidate_genes,
  up_color = group_colors[["DN_G"]],
  down_color = group_colors[["DN_T"]]
)

ggsave(
  "11_MA_DN_G_vs_Control_G.png",
  plot = p_ma_G,
  width = 8,
  height = 6,
  bg = "white"
)

ggsave(
  "12_MA_DN_T_vs_Control_T.png",
  plot = p_ma_T,
  width = 8,
  height = 6,
  bg = "white"
)

ggsave(
  "13_MA_Control_G_vs_Control_T.png",
  plot = p_ma_Control_GT,
  width = 8,
  height = 6,
  bg = "white"
)

ggsave(
  "14_MA_DN_G_vs_DN_T.png",
  plot = p_ma_DN_GT,
  width = 8,
  height = 6,
  bg = "white"
)


# ============================================================
# 23. DEG gene sets
# ============================================================

NC_DN_G <- DEG_G %>%
  transmute(
    gene,
    FC = log2FoldChange,
    padj
  )

NC_DN_G_up <- NC_DN_G %>%
  filter(
    FC > 1,
    padj < 0.01
  )

NC_DN_G_down <- NC_DN_G %>%
  filter(
    FC < -1,
    padj < 0.01
  )

NC_DN_T <- DEG_T %>%
  transmute(
    gene,
    FC = log2FoldChange,
    padj
  )

NC_DN_T_up <- NC_DN_T %>%
  filter(
    FC > 1,
    padj < 0.01
  )

NC_DN_T_down <- NC_DN_T %>%
  filter(
    FC < -1,
    padj < 0.01
  )

DN_DN <- DEG_DN_GT %>%
  transmute(
    gene,
    FC = log2FoldChange,
    padj
  )

DN_DN_up <- DN_DN %>%
  filter(
    FC > 1,
    padj < 0.01
  )

# Preserve original definition:
# "No differences in DN_G vs DN_T" = -1 < log2FC < 1
# There is no padj criterion here.
DN_DN_none <- DN_DN %>%
  filter(
    !is.na(FC),
    FC < 1,
    FC > -1
  )

DN_DN_down <- DN_DN %>%
  filter(
    FC < -1,
    padj < 0.01
  )

NC_NC <- DEG_Control_GT %>%
  transmute(
    gene,
    FC = log2FoldChange,
    padj
  )

NC_NC_up <- NC_NC %>%
  filter(
    FC > 1,
    padj < 0.01
  )

NC_NC_down <- NC_NC %>%
  filter(
    FC < -1,
    padj < 0.01
  )


# ============================================================
# 24. Loss of compartment-specific expression in DN
# ============================================================

loss_G_genes <- Reduce(
  intersect,
  list(
    NC_DN_G_down$gene,
    NC_NC_up$gene,
    DN_DN_none$gene
  )
)

loss_T_genes <- Reduce(
  intersect,
  list(
    NC_DN_T_down$gene,
    NC_NC_down$gene,
    DN_DN_none$gene
  )
)

cat(
  "Loss of glomerulus-specific expression:",
  length(loss_G_genes),
  "\n"
)

cat(
  "Loss of tubule-specific expression:",
  length(loss_T_genes),
  "\n"
)

write.csv(
  data.frame(gene = loss_G_genes),
  "Loss_G_specificity_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene = loss_T_genes),
  "Loss_T_specificity_genes.csv",
  row.names = FALSE
)


# ============================================================
# 25. Venn diagrams for compartment-specific loss
#
# DR Venn has been removed.
# ============================================================

ven_G <- list(
  "Down in DN_G vs Control_G" = NC_DN_G_down$gene,
  "UP in G of Control_G vs Control_T" = NC_NC_up$gene,
  "No differences in DN_G vs DN_T" = DN_DN_none$gene
)

p_ven_G <- ggVennDiagram(
  ven_G,
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#F4FAFE"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  "15_Venn_loss_G_specificity.png",
  plot = p_ven_G,
  width = 8,
  height = 7,
  bg = "white"
)


ven_T <- list(
  "Down in DN_T vs Control_T" = NC_DN_T_down$gene,
  "UP in T of Control_G vs Control_T" = NC_NC_down$gene,
  "No differences in DN_G vs DN_T" = DN_DN_none$gene
)

p_ven_T <- ggVennDiagram(
  ven_T,
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "#F4FAFE"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  "16_Venn_loss_T_specificity.png",
  plot = p_ven_T,
  width = 8,
  height = 7,
  bg = "white"
)


# ============================================================
# 26. DEG overlap: glomerulus vs tubule
#     Useful for Figure S3
# ============================================================

common_down_genes <- intersect(
  NC_DN_G_down$gene,
  NC_DN_T_down$gene
)

common_up_genes <- intersect(
  NC_DN_G_up$gene,
  NC_DN_T_up$gene
)

ven_down <- list(
  "DN_G vs Control_G" = NC_DN_G_down$gene,
  "DN_T vs Control_T" = NC_DN_T_down$gene
)

p_ven_down <- ggVennDiagram(
  ven_down,
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "orange"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  "17_Venn_downregulated_G_T_overlap.png",
  plot = p_ven_down,
  width = 7,
  height = 6,
  bg = "white"
)


ven_up <- list(
  "DN_G vs Control_G" = NC_DN_G_up$gene,
  "DN_T vs Control_T" = NC_DN_T_up$gene
)

p_ven_up <- ggVennDiagram(
  ven_up,
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(
    low = "#F4FAFE",
    high = "orange"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  "18_Venn_upregulated_G_T_overlap.png",
  plot = p_ven_up,
  width = 7,
  height = 6,
  bg = "white"
)


# ============================================================
# 27. Top 10 gene Boxplot
#
# KEPT because user requested boxplots to remain.
# This plot remains independent of Table S2.
#
# The original script manually overwrote the calculated top 10 with
# this gene list, so that behavior is preserved.
# ============================================================

top_10_genes <- c(
  "ANK1",
  "COL1A1",
  "DRICH1",
  "FN1",
  "GPC2",
  "H1-1",
  "LMAN1L",
  "MARCO",
  "POU6F2",
  "SPDYA"
)

top10 <- expression_long %>%
  filter(
    Gene %in% top_10_genes
  ) %>%
  mutate(
    Gene = factor(
      Gene,
      levels = top_10_genes
    )
  )

p_top10_box <- ggplot(
  top10,
  aes(
    x = Gene,
    y = Expression,
    fill = Group
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.5,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(color = Group),
    size = 1,
    position = position_jitterdodge(
      jitter.width = 0.07,
      jitter.height = 0.2
    )
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    x = "",
    y = "Log2 Expression (CPM)"
  ) +
  theme_minimal() +
  guides(
    fill = guide_legend(nrow = 1)
  ) +
  coord_cartesian(
    ylim = c(0, 13)
  ) +
  theme(
    axis.title = element_text(size = 23),
    axis.text.x = element_text(
      color = "black",
      size = 18,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      color = "black",
      size = 18
    ),
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  )

ggsave(
  "19_Top10_genes_boxplot.png",
  plot = p_top10_box,
  width = 8,
  height = 5,
  bg = "white"
)


# ============================================================
# 28. CPM pattern line plot
#     KEPT, but original overwrite bug is corrected
# ============================================================

dn_cols_cpm <- calculate_cpm(
  id_G[, dn_index, drop = FALSE]
)

nc_cols_cpm <- calculate_cpm(
  id_G[, nc_index, drop = FALSE]
)

dt_cols_cpm <- calculate_cpm(
  id_T[, dt_index, drop = FALSE]
)

nt_cols_cpm <- calculate_cpm(
  id_T[, nt_index, drop = FALSE]
)


expressionMean <- function(
    data
) {

  data_mean <- rowMeans(
    as.matrix(data)
  )

  data.frame(
    data = data_mean,
    gene = names(data_mean),
    row.names = NULL
  )
}


nc_mean <- expressionMean(
  nc_cols_cpm
)

dn_mean <- expressionMean(
  dn_cols_cpm
)

nt_mean <- expressionMean(
  nt_cols_cpm
)

dt_mean <- expressionMean(
  dt_cols_cpm
)

nc_mean$condition <- "Control_G"
dn_mean$condition <- "DN_G"
nt_mean$condition <- "Control_T"
dt_mean$condition <- "DN_T"

pp_int <- c(
  "AKR1A1",
  "ALDH2",
  "ALDH1B1"
)

nc_gene <- subset(
  nc_mean,
  gene %in% pp_int
)

dn_gene <- subset(
  dn_mean,
  gene %in% pp_int
)

nt_gene <- subset(
  nt_mean,
  gene %in% pp_int
)

dt_gene <- subset(
  dt_mean,
  gene %in% pp_int
)

# FIX:
# Original code created the G merge and then immediately overwrote
# it with the T merge. All four groups are retained here.
pattern_merge <- bind_rows(
  nc_gene,
  nt_gene,
  dn_gene,
  dt_gene
)

pattern_merge$condition <- factor(
  pattern_merge$condition,
  levels = c(
    "Control_G",
    "Control_T",
    "DN_G",
    "DN_T"
  )
)

p_pattern_line <- ggplot(
  pattern_merge,
  aes(
    x = condition,
    y = data,
    group = gene,
    color = gene
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  labs(
    x = " ",
    y = "Expression (CPM)",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(
      size = 20,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

ggsave(
  "20_Pattern_lineplot.png",
  plot = p_pattern_line,
  width = 8,
  height = 5,
  bg = "white"
)


# ============================================================
# 29. E0/E1 heatmap
#
# REMOVED as requested.
# ============================================================


# ============================================================
# 30. GO helper
# ============================================================

run_go <- function(
    genes,
    analysis_name
) {

  genes <- unique(
    genes
  )

  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]

  if (length(genes) == 0) {

    warning(
      "No genes available for GO analysis: ",
      analysis_name
    )

    return(NULL)
  }

  enrichGO(
    gene = genes,
    universe = rownames(id),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}


save_go_plots <- function(
    go_object,
    prefix,
    title
) {

  if (
    is.null(go_object) ||
      nrow(as.data.frame(go_object)) == 0
  ) {

    warning(
      "No significant GO terms for: ",
      title
    )

    return(invisible(NULL))
  }

  write.csv(
    as.data.frame(go_object),
    paste0(
      prefix,
      "_GO_result.csv"
    ),
    row.names = FALSE
  )

  p_dot <- dotplot(
    go_object,
    label_format = 80,
    showCategory = 20
  ) +
    ggtitle(
      title
    ) +
    theme(
      axis.title = element_text(size = 18),
      axis.text.x = element_text(
        color = "black",
        size = 14
      ),
      axis.text.y = element_text(
        color = "black",
        size = 11
      ),
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13)
    )

  ggsave(
    paste0(
      prefix,
      "_GO_dotplot.png"
    ),
    plot = p_dot,
    width = 9,
    height = 6,
    bg = "white"
  )

  p_bar <- barplot(
    go_object,
    showCategory = 20,
    label_format = 80
  ) +
    ggtitle(
      title
    ) +
    theme(
      axis.title = element_text(size = 18),
      axis.text.x = element_text(
        color = "black",
        size = 12
      ),
      axis.text.y = element_text(
        color = "black",
        size = 11
      ),
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13)
    )

  ggsave(
    paste0(
      prefix,
      "_GO_barplot.png"
    ),
    plot = p_bar,
    width = 9,
    height = 6,
    bg = "white"
  )
}


# ============================================================
# 31. GO analyses
#
# Main convergence analysis:
#   loss_G_genes
#   loss_T_genes
#
# Supplementary DN-associated downregulated analysis:
#   G_down
#   T_down
#   common_down
# ============================================================

go_loss_G <- run_go(
  loss_G_genes,
  "Loss of glomerulus-specific expression"
)

go_loss_T <- run_go(
  loss_T_genes,
  "Loss of tubule-specific expression"
)

G_only_down_genes <- setdiff(
  NC_DN_G_down$gene,
  NC_DN_T_down$gene
)

length(G_only_down_genes)

go_G_down <- run_go(
  NC_DN_G_down$gene,
  "DN glomerulus downregulated genes"
)

go_T_down <- run_go(
  NC_DN_T_down$gene,
  "DN tubule downregulated genes"
)

go_common_down <- run_go(
  common_down_genes,
  "Common downregulated genes"
)

save_go_plots(
  go_loss_G,
  "21_loss_G",
  "Loss of glomerulus-specific expression"
)

save_go_plots(
  go_loss_T,
  "22_loss_T",
  "Loss of tubule-specific expression"
)

save_go_plots(
  go_G_down,
  "23_DN_G_down",
  "DN glomerulus downregulated genes"
)

save_go_plots(
  go_T_down,
  "24_DN_T_down",
  "DN tubule downregulated genes"
)

save_go_plots(
  go_common_down,
  "25_common_down",
  "Common downregulated genes"
)


# ============================================================
# 32. KEGG
#
# KEPT because user requested remaining analyses to remain.
# The previous GO input depended on the DR Venn object.
# After DR removal, use the glomerulus-specific-loss gene set
# as the legacy KEGG input.
# ============================================================

loss_G_entrez <- bitr(
  loss_G_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

if (
  !is.null(loss_G_entrez) &&
    nrow(loss_G_entrez) > 0
) {

  enKegg <- enrichKEGG(
    loss_G_entrez$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    minGSSize = 10,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (
    !is.null(enKegg) &&
      nrow(as.data.frame(enKegg)) > 0
  ) {

    write.csv(
      as.data.frame(enKegg),
      "26_KEGG_loss_G_result.csv",
      row.names = FALSE
    )

    p_kegg <- dotplot(
      enKegg,
      label_format = 80,
      showCategory = 20
    )

    ggsave(
      "26_KEGG_dotplot.png",
      plot = p_kegg,
      width = 9,
      height = 6,
      bg = "white"
    )
  }
}


# ============================================================
# 33. Pathview
#
# KEPT from original analysis.
# Uses hsa04510 as in the original script.
# ============================================================

# if (
#   requireNamespace(
#     "pathview",
#     quietly = TRUE
#   ) &&
#     exists("loss_G_entrez") &&
#     !is.null(loss_G_entrez) &&
#     nrow(loss_G_entrez) > 0
# ) {
# 
#   pathview_df <- NC_DN_G %>%
#     filter(
#       gene %in% loss_G_genes
#     ) %>%
#     left_join(
#       loss_G_entrez,
#       by = c(
#         "gene" = "SYMBOL"
#       )
#     ) %>%
#     filter(
#       !is.na(ENTREZID)
#     ) %>%
#     distinct(
#       gene,
#       .keep_all = TRUE
#     )
# 
#   pathview_data <- pathview_df$FC
# 
#   names(pathview_data) <- as.character(
#     pathview_df$ENTREZID
#   )
# 
#   if (length(pathview_data) > 0) {
# 
#     pathview::pathview(
#       gene.data = pathview_data,
#       pathway.id = "hsa04510",
#       species = "hsa",
#       limit = list(
#         gene = max(
#           abs(pathview_data),
#           na.rm = TRUE
#         ),
#         cpd = 1
#       )
#     )
# 
#     if (
#       file.exists(
#         "hsa04510.pathview.png"
#       )
#     ) {
# 
#       file.rename(
#         "hsa04510.pathview.png",
#         "27_KEGG_pathview_hsa04510.png"
#       )
#     }
#   }
# }


# ============================================================
# 34. STRING PPI helper
# ============================================================

save_string_network <- function(
    proteins,
    output_file
) {

  proteins <- unique(
    proteins
  )

  tryCatch(
    {

      proteins_mapped <- rba_string_map_ids(
        ids = proteins,
        species = 9606
      )

      int_net <- rba_string_interactions_network(
        ids = proteins_mapped$stringId,
        species = 9606,
        required_score = 700
      )

      connected_proteins <- unique(
        c(
          int_net$preferredName_A,
          int_net$preferredName_B
        )
      )

      proteins_filtered <- proteins[
        proteins %in% connected_proteins
      ]

      proteins_mapped_filtered <- rba_string_map_ids(
        ids = proteins_filtered,
        species = 9606
      )

      if (
        file.exists(
          "string_network_image.png"
        )
      ) {
        file.remove(
          "string_network_image.png"
        )
      }

      rba_string_network_image(
        ids = proteins_mapped_filtered$stringId,
        image_format = "image",
        species = 9606,
        save_image = TRUE,
        required_score = 700,
        network_flavor = "confidence"
      )

      if (
        file.exists(
          "string_network_image.png"
        )
      ) {

        file.rename(
          "string_network_image.png",
          output_file
        )
      }
    },
    error = function(e) {

      warning(
        "STRING network failed for ",
        output_file,
        ": ",
        conditionMessage(e)
      )
    }
  )
}


# ============================================================
# 35. PPI networks
#
# Original G PPI is kept.
# T PPI is also generated because the manuscript contains both.
# ============================================================

save_string_network(
  NC_DN_G_down$gene,
  "28_PPI_network_DN_G_down.png"
)

save_string_network(
  NC_DN_T_down$gene,
  "29_PPI_network_DN_T_down.png"
)


# ============================================================
# 36. Outlier gene barplots
#
# FIX:
# Original code tried to use dt_cols after it had already been
# removed and after the outlier genes had been deleted from id.
# Use id_unfiltered here so the excluded genes remain available.
# ============================================================

outlier_cpm_unfiltered <- calculate_cpm(
  id_unfiltered
)

control_G_cols <- metaData$id[
  metaData$sample == "Control_G"
]

DN_T_cols <- metaData$id[
  metaData$sample == "DN_T"
]

out_control_G <- outlier_cpm_unfiltered[
  intersect(
    c(
      "FOXN2",
      "GPX3"
    ),
    rownames(outlier_cpm_unfiltered)
  ),
  control_G_cols,
  drop = FALSE
]

out_control_G_long <- as.data.frame(
  out_control_G
) %>%
  rownames_to_column(
    "Gene"
  ) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  )

p_outlier_G <- ggplot(
  out_control_G_long,
  aes(
    x = Sample,
    y = Expression,
    fill = Gene
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "",
    y = "Expression (CPM)",
    x = " "
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_text(
      size = 15,
      angle = 30,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 14,
      colour = "black"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )


out_DN_T <- outlier_cpm_unfiltered[
  intersect(
    c(
      "IGHG1",
      "IGHG2",
      "IGHG3",
      "IGHG4",
      "IGKC",
      "IGLL5"
    ),
    rownames(outlier_cpm_unfiltered)
  ),
  DN_T_cols,
  drop = FALSE
]

out_DN_T_long <- as.data.frame(
  out_DN_T
) %>%
  rownames_to_column(
    "Gene"
  ) %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  )

p_outlier_T <- ggplot(
  out_DN_T_long,
  aes(
    x = Sample,
    y = Expression,
    fill = Gene
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    title = "",
    y = "Expression (CPM)",
    x = " "
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    axis.text.x = element_text(
      size = 15,
      angle = 30,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 14,
      colour = "black"
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

p_outlier_combined <- ggpubr::ggarrange(
  p_outlier_G,
  p_outlier_T,
  ncol = 1,
  nrow = 2,
  labels = c(
    "A",
    "B"
  )
)

ggsave(
  "30_Barplot_outlier_genes.png",
  plot = p_outlier_combined,
  width = 11,
  height = 10,
  bg = "white"
)

ggsave(
  "30A_Barplot_outlier_Control_G.png",
  plot = p_outlier_G,
  width = 11,
  height = 5,
  bg = "white"
)

ggsave(
  "30B_Barplot_outlier_DN_T.png",
  plot = p_outlier_T,
  width = 11,
  height = 5,
  bg = "white"
)
# ============================================================
# 37. Euler diagrams using ACTUAL gene sets
#
# The original code manually entered counts.
# Use the actual gene sets so the result updates automatically.
# ============================================================

fit_G <- euler(
  list(
    "Down in DN_G vs Control_G" = NC_DN_G_down$gene,
    "UP in G of Control_G vs Control_T" = NC_NC_up$gene,
    "No differences in DN_G vs DN_T" = DN_DN_none$gene
  )
)

png(
  "31_Euler_diagram_loss_G.png",
  width = 900,
  height = 700
)

plot(
  fit_G,
  quantities = TRUE,
  labels = list(
    font = 2
  )
)

dev.off()


fit_T <- euler(
  list(
    "Down in DN_T vs Control_T" = NC_DN_T_down$gene,
    "UP in T of Control_G vs Control_T" = NC_NC_down$gene,
    "No differences in DN_G vs DN_T" = DN_DN_none$gene
  )
)

png(
  "32_Euler_diagram_loss_T.png",
  width = 900,
  height = 700
)

plot(
  fit_T,
  quantities = TRUE,
  labels = list(
    font = 2
  )
)

dev.off()


# ============================================================
# 38. Reproduction checks
# ============================================================

reproduction_check <- data.frame(
  Result = c(
    "DN_G down",
    "DN_G up",
    "DN_T down",
    "DN_T up",
    "Common down",
    "Common up",
    "Control G-high",
    "Control T-high",
    "DN G-high",
    "DN T-high",
    "Loss G specificity",
    "Loss T specificity"
  ),
  Observed = c(
    nrow(NC_DN_G_down),
    nrow(NC_DN_G_up),
    nrow(NC_DN_T_down),
    nrow(NC_DN_T_up),
    length(common_down_genes),
    length(common_up_genes),
    nrow(NC_NC_up),
    nrow(NC_NC_down),
    nrow(DN_DN_up),
    nrow(DN_DN_down),
    length(loss_G_genes),
    length(loss_T_genes)
  ),
  Thesis_reported = c(
    625,
    97,
    196,
    3,
    74,
    0,
    390,
    171,
    53,
    21,
    168,
    67
  )
)

reproduction_check$Match <- with(
  reproduction_check,
  Observed == Thesis_reported
)

print(
  reproduction_check
)

write.csv(
  reproduction_check,
  "Reproduction_check.csv",
  row.names = FALSE
)


# ============================================================
# 39. Save useful R objects and session information
# ============================================================

saveRDS(
  dds,
  "Study2_dds.rds"
)

saveRDS(
  list(
    DEG_G = DEG_G,
    DEG_T = DEG_T,
    DEG_Control_GT = DEG_Control_GT,
    DEG_DN_GT = DEG_DN_GT,
    loss_G_genes = loss_G_genes,
    loss_T_genes = loss_T_genes,
    common_down_genes = common_down_genes,
    common_up_genes = common_up_genes
  ),
  "Study2_analysis_results.rds"
)

capture.output(
  sessionInfo(),
  file = "sessionInfo.txt"
)

cat(
  "\nAnalysis complete.\n"
)
