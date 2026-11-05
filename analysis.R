#!/usr/bin/env Rscript
# Spatial transcriptomics of human diabetic nephropathy: revised manuscript.
# Entry point: Rscript analysis.R --help
# Unnormalized gene counts -> AOI voom-dream -> patient-compartment pseudobulk
# -> interaction/attenuation -> selected-gene RMS -> GO -> STRING -> figures.
# All 46 AOIs and all nonzero-total genes are retained in the primary analysis.

options(stringsAsFactors = FALSE, timeout = max(600, getOption("timeout")))

parse_options <- function(args) {
  opt <- list(download = FALSE, check_only = FALSE, skip_go = FALSE,
              skip_ppi = FALSE, skip_plots = FALSE, sensitivity = FALSE,
              string_live = FALSE, workers = 1L, output = NULL,
              counts = NULL, metadata = NULL, help = FALSE)
  flags <- c("download", "check-only", "skip-go", "skip-ppi", "skip-plots",
             "sensitivity", "string-live", "help")
  for (a in args) {
    key <- sub("=.*$", "", sub("^--", "", a))
    if (!startsWith(a, "--")) stop("Unknown argument: ", a)
    if (key %in% flags && !grepl("=", a, fixed = TRUE)) {
      opt[[gsub("-", "_", key, fixed = TRUE)]] <- TRUE
    } else if (key %in% c("workers", "output", "counts", "metadata") &&
               grepl("=", a, fixed = TRUE)) {
      opt[[key]] <- sub("^[^=]*=", "", a)
    } else stop("Unknown argument: ", a, ". Use --help.")
  }
  n <- suppressWarnings(as.numeric(opt$workers))
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n))
    stop("--workers must be a positive integer.")
  opt$workers <- as.integer(n)
  opt
}

need_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                            ". See README.md; no packages are installed automatically.")
}

record_versions <- function(pkgs, reference_dir, out) {
  tracked <- unique(c("readxl", pkgs, "lme4", "lmerTest", "Matrix", "minqa", "nloptr",
                      "Rcpp", "RhpcBLASctl", "GO.db", "AnnotationDbi", "org.Hs.eg.db"))
  tracked <- intersect(tracked, rownames(utils::installed.packages()))
  pv <- data.frame(Package = c("R", tracked), Version = c(as.character(getRversion()),
    vapply(tracked, function(p) as.character(utils::packageVersion(p)), character(1))))
  original <- utils::read.csv(file.path(reference_dir, "original_package_versions.csv"))
  # The original session report also records model/optimizer dependencies whose
  # versions can affect whether a borderline model converges successfully.
  si <- paste(readLines(file.path(reference_dir, "original_sessionInfo.txt"), warn = FALSE), collapse = " ")
  tokens <- regmatches(si, gregexpr("[A-Za-z][A-Za-z0-9.]*_[0-9][0-9.-]*", si))[[1]]
  deps <- data.frame(Package = sub("_.*$", "", tokens), Version = sub("^[^_]+_", "", tokens))
  original <- rbind(original, deps, data.frame(Package = "R", Version = "4.3.0"))
  original <- original[!duplicated(original$Package), , drop = FALSE]
  original <- original[original$Package %in% pv$Package, , drop = FALSE]
  original$Version <- vapply(original$Version, function(v) as.character(package_version(v)), character(1))
  vc <- merge(pv, original, by = "Package", suffixes = c("_run", "_original"))
  vc$Matches <- vc$Version_run == vc$Version_original
  write_table(pv, file.path(out, "package_versions.csv"))
  write_table(vc, file.path(out, "package_version_comparison.csv"))
  if (any(!vc$Matches)) warning("Package versions differ from the original run; review package_version_comparison.csv.")
  invisible(vc)
}

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "NA")
}

write_matrix <- function(x, path, id = "gene") {
  out <- data.frame(ID = rownames(x), x, check.names = FALSE, row.names = NULL)
  names(out)[1] <- id
  write_table(out, path)
}

canonical_aoi <- function(x) {
  # Accepted source examples: NC | Normal-1-G1, DN | DN-2-G1,
  # Normal-1-G1, Control-1-G1, DN-2-G1. Never match by column position.
  x <- trimws(as.character(x))
  x <- sub("^(NC|DN)[[:space:]]*\\|[[:space:]]*", "", x)
  x <- sub("^Normal-", "Control-", x)
  if (anyNA(x) || any(!grepl("^(Control|DN)-[1-4]-[GT][1-9][0-9]*$", x)))
    stop("Unrecognized AOI label(s): ", paste(x[is.na(x) |
           !grepl("^(Control|DN)-[1-4]-[GT][1-9][0-9]*$", x)], collapse = ", "))
  x
}

read_inputs <- function(count_file, metadata_file) {
  need_packages("readxl")
  con <- if (grepl("\\.gz$", count_file, ignore.case = TRUE)) gzfile(count_file, open = "rt") else file(count_file, open = "rt")
  on.exit(close(con), add = TRUE)
  d <- utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) < 2L || nrow(d) == 0L) stop("Empty or invalid raw count matrix.")
  genes <- as.character(d[[1]])
  if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) || any(trimws(genes) != genes))
    stop("Gene symbols must be nonempty, unique, and free of surrounding whitespace.")
  if (!all(vapply(d[-1], is.numeric, logical(1))))
    stop("All columns after the first gene-symbol column must be numeric counts.")
  counts <- as.matrix(d[-1]); rownames(counts) <- genes
  if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8))
    stop("Counts must be finite, nonnegative integers. Do not use logCPM or normalized values.")
  old_cols <- colnames(counts)
  colnames(counts) <- canonical_aoi(old_cols)
  if (anyDuplicated(colnames(counts))) stop("Duplicate count AOIs after name standardization.")
  meta <- as.data.frame(readxl::read_excel(metadata_file, sheet = 1), stringsAsFactors = FALSE)
  required <- c("SlideName", "ScanLabel", "ROILabel", "AOISurfaceArea",
                "AOINucleiCount", "ROICoordinateX", "ROICoordinateY")
  if (!all(required %in% names(meta))) stop("Metadata must contain: ", paste(required, collapse = ", "))
  meta$AOI_ID <- canonical_aoi(meta$ROILabel)
  if (anyDuplicated(meta$AOI_ID)) stop("Duplicate AOIs in slide assignment.")
  if (!setequal(colnames(counts), meta$AOI_ID)) {
    stop("AOI mismatch. Counts only: ", paste(setdiff(colnames(counts), meta$AOI_ID), collapse = ", "),
         "; metadata only: ", paste(setdiff(meta$AOI_ID, colnames(counts)), collapse = ", "))
  }
  meta$Disease <- sub("-.*$", "", meta$AOI_ID)
  meta$Patient_ID <- sub("-[GT][0-9]+$", "", meta$AOI_ID)
  meta$Patient_Number <- as.integer(sub("^(Control|DN)-([0-9]+)-.*$", "\\2", meta$AOI_ID))
  meta$Compartment <- sub("^.*-([GT])[0-9]+$", "\\1", meta$AOI_ID)
  meta$AOI_Number <- as.integer(sub("^.*-[GT]([0-9]+)$", "\\1", meta$AOI_ID))
  meta$Group <- paste(meta$Disease, meta$Compartment, sep = "_")
  meta$Patient_Compartment <- paste(meta$Patient_ID, meta$Compartment, sep = "__")
  scan_disease <- ifelse(meta$ScanLabel %in% c("NC", "Normal", "Control"), "Control",
                         ifelse(meta$ScanLabel == "DN", "DN", NA_character_))
  if (anyNA(scan_disease) || any(scan_disease != meta$Disease))
    stop("ScanLabel and ROILabel disagree on disease group.")
  if (anyNA(meta$SlideName) || any(!nzchar(trimws(meta$SlideName)))) stop("Missing slide assignment.")
  for (nm in c("AOISurfaceArea", "AOINucleiCount", "ROICoordinateX", "ROICoordinateY")) {
    if (!is.numeric(meta[[nm]]) || any(!is.finite(meta[[nm]]))) stop("Invalid metadata field: ", nm)
  }
  if (any(meta$AOISurfaceArea <= 0) || any(meta$AOINucleiCount < 0)) stop("Invalid AOI area/nuclei counts.")
  # Stable order for numerical reproducibility, independent of input order.
  ii <- order(match(meta$Disease, c("Control", "DN")), meta$Patient_Number,
              match(meta$Compartment, c("G", "T")), meta$AOI_Number)
  meta <- meta[ii, , drop = FALSE]
  rownames(meta) <- meta$AOI_ID
  counts <- counts[, meta$AOI_ID, drop = FALSE]
  meta$Original_count_column <- old_cols[match(meta$AOI_ID, canonical_aoi(old_cols))]
  meta$Disease <- factor(meta$Disease, levels = c("Control", "DN"))
  meta$Compartment <- factor(meta$Compartment, levels = c("T", "G"))
  meta$Group <- factor(meta$Group, levels = c("Control_T", "Control_G", "DN_T", "DN_G"))
  meta$Patient_ID <- factor(meta$Patient_ID)
  meta$Patient_Compartment <- factor(meta$Patient_Compartment)
  expected <- c(Control_G = 12L, Control_T = 12L, DN_G = 10L, DN_T = 12L)
  actual <- table(meta$Group)[names(expected)]
  if (nrow(meta) != 46L || !identical(as.integer(actual), unname(expected)) ||
      length(unique(meta$Patient_ID)) != 8L) stop("This release expects the 46-AOI, eight-patient manuscript cohort.")
  pair <- table(meta$Patient_ID, meta$Compartment)
  incomplete <- rownames(pair)[rowSums(pair > 0) < 2L]
  if (!identical(incomplete, "DN-1") || pair["DN-1", "G"] != 0L)
    stop("Unexpected patient-compartment coverage; DN-1 should be the only missing G profile.")
  meta$Low_nuclei_flag <- meta$AOINucleiCount < 200
  meta$Low_area_flag <- meta$AOISurfaceArea < 16000
  meta$Included_in_primary <- TRUE
  # No read/probe QC flags are inferred from counts or invented from this workbook.
  filter <- data.frame(gene = genes, Total_raw_count = rowSums(counts), Kept = rowSums(counts) > 0)
  counts <- counts[filter$Kept, , drop = FALSE]
  if (any(colSums(counts) <= 0)) stop("An AOI has zero library size.")
  list(counts = counts, meta = meta, filter = filter)
}

fetch_geo <- function(data_dir) {
  need_packages("GEOquery")
  message("Downloading the two supplementary input files from GSE346609, if publicly available.")
  # GEO may prepend the accession and/or gzip the CSV. Discover filenames rather
  # than assuming an unverified direct supplementary-file URL.
  GEOquery::getGEOSuppFiles("GSE346609", makeDirectory = FALSE, baseDir = data_dir,
    filter_regex = "(raw_counts\\.csv(\\.gz)?|Slide_Assignment\\.xlsx)$")
}

resolve_input <- function(explicit, data_dir, pattern, default) {
  if (!is.null(explicit)) {
    if (!file.exists(explicit)) stop("Input does not exist: ", explicit)
    return(normalizePath(explicit))
  }
  candidates <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  preferred <- file.path(data_dir, default)
  if (file.exists(preferred)) {
    if (length(candidates) > 1L) warning("Multiple input versions found; using ", preferred,
      ". Use --counts= or --metadata= to select a different file.")
    return(preferred)
  }
  if (length(candidates) != 1L) stop("Place ", default, " in data/ (or use --download after GEO release). ",
                                    "If multiple versions exist, specify --counts= or --metadata= explicitly.")
  candidates
}

contrast_definitions <- c(
  G_DN_vs_Control = "GroupDN_G - GroupControl_G",
  T_DN_vs_Control = "GroupDN_T - GroupControl_T",
  Control_GT = "GroupControl_G - GroupControl_T",
  DN_GT = "GroupDN_G - GroupDN_T",
  Interaction = "(GroupDN_G - GroupDN_T) - (GroupControl_G - GroupControl_T)")

fit_dream <- function(counts, meta, form, out, workers = 1L, universe = rownames(counts)) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  stopifnot(identical(colnames(counts), rownames(meta)))
  meta <- droplevels(meta)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  if (!nrow(counts) || any(colSums(counts) == 0)) stop("Empty model input.")
  design <- model.matrix(~0 + Group, meta)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient fixed design.")
  if (any(table(unique(meta[, c("Patient_ID", "Group")])$Group) < 2L))
    stop("Fewer than two biological replicates in a group.")
  bp <- if (workers == 1L) BiocParallel::SerialParam(progressbar = FALSE) else
    BiocParallel::SnowParam(workers, type = "SOCK", progressbar = FALSE)
  if (workers > 1L) on.exit(BiocParallel::bpstop(bp), add = TRUE)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
  L <- variancePartition::makeContrastsDream(form, meta, contrasts = contrast_definitions)
  write_matrix(L, file.path(out, "contrast_matrix.csv"), "coefficient")
  write_table(meta, file.path(out, "metadata.csv"))
  write_table(data.frame(Profile_ID = colnames(counts), dge$samples, row.names = NULL),
              file.path(out, "TMM_library_factors.csv"))
  message("Estimating voom weights for ", nrow(counts), " genes and ", ncol(counts), " profiles.")
  v <- variancePartition::voomWithDreamWeights(dge, form, meta,
    normalize.method = "none", prior.count = 0.5, rescaleWeightsAfter = FALSE,
    plot = FALSE, save.plot = TRUE, BPPARAM = bp)
  message("Fitting dream contrasts and gene-specific degrees of freedom.")
  f0 <- variancePartition::dream(v, form, meta, L = L, ddf = "Satterthwaite",
                                  REML = TRUE, BPPARAM = bp)
  utils::capture.output(list(Formula = form, Voom_errors = attr(v, "errors"),
    Dream_errors = attr(f0, "errors"), Initial_errors = attr(f0, "error.initial")),
    file = file.path(out, "fit_diagnostics.txt"))
  fit <- variancePartition::eBayes(f0)
  ans <- list()
  for (nm in names(contrast_definitions)) {
    # Directly use the gene/contrast-specific moderated df. In 1.32.5, the
    # topTable(confint=TRUE) route is unsuitable for these multiple contrasts.
    tt <- variancePartition::topTable(fit, coef = nm, number = Inf,
      sort.by = "none", adjust.method = "BH", p.value = 1, lfc = 0, confint = FALSE)
    gi <- match(rownames(tt), rownames(fit$coefficients))
    ci <- match(nm, colnames(fit$coefficients))
    stopifnot(!anyNA(gi), !is.na(ci))
    s2 <- if (is.matrix(fit$s2.post)) fit$s2.post[gi, ci] else fit$s2.post[gi]
    df <- if (is.matrix(fit$df.total)) fit$df.total[gi, ci] else fit$df.total[gi]
    se <- fit$stdev.unscaled[gi, ci] * sqrt(s2)
    margin <- stats::qt(0.975, df) * se
    compact <- data.frame(gene = rownames(tt), log2FoldChange = tt$logFC,
      lfcSE = se, df = df, stat = tt$t, pvalue = tt$P.Value, padj = tt$adj.P.Val,
      CI_low = tt$logFC - margin, CI_high = tt$logFC + margin, AveExpr = tt$AveExpr)
    full <- compact[match(universe, compact$gene), , drop = FALSE]
    full$gene <- universe; rownames(full) <- NULL
    full$Present_in_run <- universe %in% rownames(counts)
    full$Fit_available <- is.finite(full$log2FoldChange) & is.finite(full$pvalue)
    full$Significant <- is.finite(full$padj) & full$padj < 0.05 & is.finite(full$log2FoldChange)
    if (nm != "Interaction") full$Significant <- full$Significant & abs(full$log2FoldChange) > 1
    full$Comparison <- nm
    full$Direction <- "Not_significant"
    full$Direction[full$Significant & full$log2FoldChange > 0] <- if (nm %in% c("Control_GT", "DN_GT")) "G_high" else "Up_in_DN"
    full$Direction[full$Significant & full$log2FoldChange < 0] <- if (nm %in% c("Control_GT", "DN_GT")) "T_high" else "Down_in_DN"
    if (nm == "Interaction") full$Direction[full$Significant] <- "Interaction"
    full$CI_type <- "95% pointwise moderated-t; not multiplicity-adjusted"
    ans[[nm]] <- full
    write_table(full, file.path(out, paste0(nm, "_ALL.csv")))
    write_table(full[full$Significant, , drop = FALSE], file.path(out, paste0(nm, "_SIGNIFICANT.csv")))
  }
  write_matrix(v$E, file.path(out, "TMM_logCPM.csv"))
  list(results = ans, logCPM = v$E)
}

classify_interaction <- function(res) {
  x <- res$Interaction
  for (g in c("Control", "DN")) {
    y <- res[[paste0(g, "_GT")]]; ii <- match(x$gene, y$gene)
    x[[paste0(g, "_GT_LFC")]] <- y$log2FoldChange[ii]
    x[[paste0(g, "_GT_FDR")]] <- y$padj[ii]
  }
  C <- x$Control_GT_LFC; D <- x$DN_GT_LFC
  sig <- is.finite(x$padj) & x$padj < 0.05 & is.finite(x$log2FoldChange)
  baseline <- is.finite(x$Control_GT_FDR) & x$Control_GT_FDR < 0.05 & is.finite(C)
  x$Class <- "Not_significant"
  x$Class[sig] <- "Other_interaction"
  x$Class[which(sig & baseline & C > 1 & is.finite(D) & D >= 0 & D < C)] <- "G_attenuation"
  x$Class[which(sig & baseline & C < -1 & is.finite(D) & D <= 0 & D > C)] <- "T_attenuation"
  # Reversals remain an annotation within Other_interaction, not a fifth class.
  x$Estimated_sign_reversal <- is.finite(C) & is.finite(D) & C * D < 0
  x$Small_control_contrast <- is.finite(C) & abs(C) < 1
  x
}

aggregate_profiles <- function(counts, expr, meta) {
  ids <- unique(as.character(meta$Patient_Compartment))
  m <- meta[match(ids, as.character(meta$Patient_Compartment)), , drop = FALSE]
  m$Profile_ID <- ids; m$AOI_ID <- NULL; m$Original_count_column <- NULL
  rownames(m) <- ids
  summed <- means <- matrix(NA_real_, nrow(counts), length(ids), dimnames = list(rownames(counts), ids))
  n_aoi <- integer(length(ids))
  for (i in seq_along(ids)) {
    keep <- as.character(meta$Patient_Compartment) == ids[i]
    summed[, i] <- rowSums(counts[, keep, drop = FALSE])
    means[, i] <- rowMeans(expr[, keep, drop = FALSE])
    n_aoi[i] <- sum(keep)
  }
  # Pseudobulk metadata should not imply that a single AOI's coordinates,
  # area, nuclei count or source label describes the aggregated profile.
  m <- m[, c("Profile_ID", "Patient_ID", "Disease", "Compartment", "Group", "Patient_Compartment", "SlideName")]
  m$n_AOI <- n_aoi
  stopifnot(isTRUE(all.equal(rowSums(summed), rowSums(counts), check.attributes = FALSE)))
  list(counts = summed, means = means, meta = m)
}

selected_distance <- function(means, meta, res, out) {
  genes <- res$Control_GT$gene[res$Control_GT$Significant]
  pairs <- table(meta$Patient_ID, meta$Compartment)
  complete <- rownames(pairs)[rowSums(pairs > 0) == 2L]
  keep <- as.character(meta$Patient_ID) %in% complete
  m <- meta[keep, , drop = FALSE]
  if (!length(genes) || !all(genes %in% rownames(means))) stop("Invalid selected-gene RMS input.")
  e <- means[genes, rownames(m), drop = FALSE]
  center <- rowMeans(e); s <- apply(e, 1, stats::sd)
  if (any(!is.finite(e)) || any(!is.finite(s) | s <= 0)) stop("Selected genes have unavailable/constant profiles.")
  z <- sweep(sweep(e, 1, center, "-"), 1, s, "/")
  d <- data.frame(Patient_ID = complete, Disease = NA_character_, RMS_zlogCPM = NA_real_, n_genes = length(genes))
  for (i in seq_along(complete)) {
    ii <- as.character(m$Patient_ID) == complete[i]
    g <- rownames(m)[ii & m$Compartment == "G"]
    t <- rownames(m)[ii & m$Compartment == "T"]
    stopifnot(length(g) == 1L, length(t) == 1L)
    d$Disease[i] <- as.character(unique(m$Disease[ii]))
    d$RMS_zlogCPM[i] <- sqrt(mean((z[, g] - z[, t])^2))
  }
  sm <- do.call(rbind, lapply(c("Control", "DN"), function(g) {
    v <- d$RMS_zlogCPM[d$Disease == g]
    data.frame(Disease = g, n_patients = length(v), Mean_RMS = mean(v), SD_RMS = stats::sd(v))
  }))
  write_table(d, file.path(out, "patient_RMS_distances.csv"))
  write_table(sm, file.path(out, "RMS_descriptive_summary.csv"))
  write_table(data.frame(gene = genes, mean = center, sd = s), file.path(out, "RMS_gene_standardization.csv"))
  write_matrix(z, file.path(out, "RMS_complete_pair_gene_zscores.csv"))
  write_table(data.frame(Patient_ID = rownames(pairs), Included = rownames(pairs) %in% complete),
              file.path(out, "RMS_patient_inclusion.csv"))
  # No inferential P value: selection and distance assessment share the cohort.
  list(patient = d, summary = sm)
}

marker_definitions <- data.frame(
  gene = c("NPHS1", "NPHS2", "PODXL", "SYNPO", "LRP2", "SLC34A1", "SLC9A3",
           "SLC12A3", "AQP2", "COL1A1", "COL3A1", "DCN", "LUM", "PTPRC", "CD68", "LST1", "CD3D"),
  Marker_group = rep(c("Podocyte-associated", "Proximal tubule", "Other tubular segments",
                       "ECM/interstitial", "Immune-associated"), c(4, 3, 2, 4, 4)))

prepare_markers <- function(expr, meta, means, pm, out) {
  genes <- marker_definitions$gene
  if (!all(genes %in% rownames(expr))) stop("Missing marker genes: ", paste(setdiff(genes, rownames(expr)), collapse = ", "))
  ref <- means[genes, , drop = FALSE]
  mu <- rowMeans(ref); sd <- apply(ref, 1, stats::sd)
  if (any(!is.finite(sd) | sd <= 0)) stop("Invalid marker reference SD.")
  a <- expr[genes, , drop = FALSE]
  za <- sweep(sweep(a, 1, mu, "-"), 1, sd, "/")
  zm <- sweep(sweep(ref, 1, mu, "-"), 1, sd, "/")
  write_table(data.frame(gene = genes, Mean_patient_compartment_logCPM = mu,
                         SD_patient_compartment_logCPM = sd), file.path(out, "marker_scaling_reference.csv"))
  write_matrix(a, file.path(out, "marker_AOI_logCPM.csv"))
  write_matrix(ref, file.path(out, "marker_patient_compartment_mean_logCPM.csv"))
  write_matrix(za, file.path(out, "marker_AOI_zscores.csv"))
  write_matrix(zm, file.path(out, "marker_patient_compartment_zscores.csv"))
  list(aoi = a, mean = ref, aoi_z = za, mean_z = zm)
}

run_go <- function(res, interaction, out, reference_dir) {
  need_packages(c("clusterProfiler", "AnnotationDbi", "org.Hs.eg.db"))
  db <- org.Hs.eg.db::org.Hs.eg.db
  genes <- unique(unlist(lapply(res, function(t) t$gene[is.finite(t$pvalue)])))
  raw <- AnnotationDbi::select(db, keys = intersect(genes, AnnotationDbi::keys(db, keytype = "SYMBOL")),
                              columns = "ENTREZID", keytype = "SYMBOL")
  raw <- unique(raw[, c("SYMBOL", "ENTREZID")])
  valid <- raw[!is.na(raw$ENTREZID), , drop = FALSE]
  ns <- table(valid$SYMBOL); ne <- table(valid$ENTREZID)
  accepted <- valid[ns[valid$SYMBOL] == 1 & ne[valid$ENTREZID] == 1, , drop = FALSE]
  write_table(raw, file.path(out, "symbol_to_entrez_all.csv"))
  write_table(accepted, file.path(out, "symbol_to_entrez_accepted.csv"))
  write_table(AnnotationDbi::metadata(db), file.path(out, "OrgDb_metadata.csv"))
  foreground <- list(
    G_attenuation = interaction$gene[interaction$Class == "G_attenuation"],
    T_attenuation = interaction$gene[interaction$Class == "T_attenuation"],
    G_Down_in_DN = res$G_DN_vs_Control$gene[res$G_DN_vs_Control$Direction == "Down_in_DN"],
    T_Down_in_DN = res$T_DN_vs_Control$gene[res$T_DN_vs_Control$Direction == "Down_in_DN"])
  backgrounds <- list(
    G_attenuation = interaction$gene[is.finite(interaction$pvalue)],
    T_attenuation = interaction$gene[is.finite(interaction$pvalue)],
    G_Down_in_DN = res$G_DN_vs_Control$gene[is.finite(res$G_DN_vs_Control$pvalue)],
    T_Down_in_DN = res$T_DN_vs_Control$gene[is.finite(res$T_DN_vs_Control$pvalue)])
  display_ids <- utils::read.csv(file.path(reference_dir, "GO_displayed_terms.csv"), check.names = FALSE)
  results <- displayed <- list(); audits <- list()
  for (nm in names(foreground)) {
    fg <- unique(accepted$ENTREZID[accepted$SYMBOL %in% foreground[[nm]]])
    bg <- unique(accepted$ENTREZID[accepted$SYMBOL %in% backgrounds[[nm]]])
    stopifnot(all(fg %in% bg))
    write_table(data.frame(gene = foreground[[nm]]), file.path(out, paste0(nm, "_input_genes.csv")))
    eg <- if (length(fg) >= 3L && length(bg)) clusterProfiler::enrichGO(
      gene = fg, universe = bg, OrgDb = db, keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 10, maxGSSize = 500, readable = TRUE) else NULL
    tab <- if (is.null(eg)) data.frame(ID = character(), Description = character(),
      GeneRatio = character(), BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
      qvalue = numeric(), geneID = character(), Count = integer()) else as.data.frame(methods::slot(eg, "result"))
    effective_bg <- if (is.null(eg)) character() else as.character(methods::slot(eg, "universe"))
    effective_fg <- intersect(fg, effective_bg)
    sig <- tab[which(is.finite(tab$p.adjust) & tab$p.adjust < 0.05), , drop = FALSE]
    desired <- display_ids$ID[display_ids$Gene_set == nm]
    show <- sig[match(desired[desired %in% sig$ID], sig$ID), , drop = FALSE]
    missing <- setdiff(desired, sig$ID)
    if (length(missing)) warning(nm, ": archived display terms absent/non-significant in this run: ", paste(missing, collapse = ", "))
    audits[[nm]] <- data.frame(Gene_set = nm, Input_n = length(foreground[[nm]]),
      Mapped_foreground_n = length(fg), Effective_foreground_n = length(effective_fg),
      Effective_background_n = length(effective_bg), Tested_terms = nrow(tab),
      Significant_terms = nrow(sig), Displayed_terms = nrow(show),
      Missing_archived_display_terms = paste(missing, collapse = ";"))
    for (what in c("ALL", "SIGNIFICANT", "DISPLAYED")) {
      x <- switch(what, ALL = tab, SIGNIFICANT = sig, DISPLAYED = show)
      write_table(x, file.path(out, paste0(nm, "_GO_BP_", what, ".csv")))
    }
    write_table(data.frame(ENTREZID = effective_bg), file.path(out, paste0(nm, "_effective_background.csv")))
    write_table(data.frame(ENTREZID = effective_fg), file.path(out, paste0(nm, "_effective_foreground.csv")))
    results[[nm]] <- tab; displayed[[nm]] <- show
  }
  audit <- do.call(rbind, audits); rownames(audit) <- NULL
  write_table(audit, file.path(out, "GO_audit.csv"))
  list(all = results, displayed = displayed, audit = audit)
}

run_ppi <- function(interaction, out, reference_dir, live = FALSE) {
  if (live) need_packages("rbioapi")
  networks <- list(); audits <- list()
  for (comp in c("G", "T")) {
    cls <- paste0(comp, "_attenuation")
    genes <- interaction$gene[interaction$Class == cls]
    if (live) {
      mapping <- rbioapi::rba_string_map_ids(ids = genes, species = 9606,
        echo_query = TRUE, timeout = 120, skip_error = FALSE, save_file = FALSE)
    } else {
      mapping <- utils::read.csv(file.path(reference_dir, paste0("PPI_", cls, "_mapping.csv")), check.names = FALSE)
      missing <- setdiff(genes, mapping$queryItem)
      if (length(missing)) stop("Archived STRING cache does not cover new ", cls,
        " genes: ", paste(missing, collapse = ", "), ". Inspect model differences or use --string-live.")
      mapping <- mapping[mapping$queryItem %in% genes, , drop = FALSE]
    }
    required <- c("queryItem", "preferredName", "stringId", "ncbiTaxonId")
    if (!all(required %in% names(mapping))) stop("Unexpected STRING mapping schema.")
    mapping$Accepted <- !is.na(mapping$stringId) & mapping$ncbiTaxonId == 9606 &
      (mapping$queryItem == mapping$preferredName |
       (mapping$queryItem == "G6PC1" & mapping$preferredName == "G6PC"))
    mapping$Accepted[is.na(mapping$Accepted)] <- FALSE
    mapping$Reason <- ifelse(mapping$Accepted, "Exact_symbol_or_verified_G6PC1_alias", "Rejected_identifier_or_species_mismatch")
    accepted <- mapping[mapping$Accepted, , drop = FALSE]
    # VEGFA -> COL18A1 is rejected by the symbol rule, not relabelled as VEGFA.
    if (anyDuplicated(accepted$queryItem) || anyDuplicated(accepted$stringId))
      stop("Ambiguous accepted STRING mapping; inspect mapping audit.")
    ids <- accepted$stringId
    if (live) {
      edges <- if (length(ids) >= 2L) rbioapi::rba_string_interactions_network(
        ids = ids, species = 9606, required_score = 700, add_nodes = 0,
        network_type = "functional", timeout = 120, skip_error = FALSE, save_file = FALSE) else data.frame()
    } else edges <- utils::read.csv(file.path(reference_dir, paste0("PPI_", cls, "_edges.csv")), check.names = FALSE)
    if (!nrow(edges)) edges <- data.frame(stringId_A = character(), stringId_B = character(), score = numeric())
    if (!all(c("stringId_A", "stringId_B", "score") %in% names(edges))) stop("Unexpected STRING edge schema.")
    edges <- edges[edges$stringId_A %in% ids & edges$stringId_B %in% ids &
      is.finite(edges$score) & edges$score >= 0.7 & edges$stringId_A != edges$stringId_B, , drop = FALSE]
    key <- paste(pmin(edges$stringId_A, edges$stringId_B), pmax(edges$stringId_A, edges$stringId_B), sep = "|")
    edges <- edges[!duplicated(key), , drop = FALSE]
    edges$gene_A <- accepted$queryItem[match(edges$stringId_A, accepted$stringId)]
    edges$gene_B <- accepted$queryItem[match(edges$stringId_B, accepted$stringId)]
    connected <- unique(c(edges$gene_A, edges$gene_B))
    write_table(mapping, file.path(out, paste0(cls, "_mapping_audit.csv")))
    write_table(edges, file.path(out, paste0(cls, "_edges.csv")))
    write_table(data.frame(gene = genes, Mapped = genes %in% accepted$queryItem,
      Connected = genes %in% connected), file.path(out, paste0(cls, "_nodes.csv")))
    audits[[cls]] <- data.frame(Gene_set = cls, Input_n = length(genes), Mapped_n = nrow(accepted),
      Connected_n = length(connected), Edges_n = nrow(edges),
      Source = if (live) "Live_STRING_query" else "Archived_STRING_mapping_and_edges",
      Confidence = 0.7, Added_nodes = 0)
    networks[[cls]] <- edges
  }
  audit <- do.call(rbind, audits); rownames(audit) <- NULL
  write_table(audit, file.path(out, "STRING_audit.csv"))
  if (live) utils::capture.output(rbioapi::rba_string_version(), file = file.path(out, "STRING_version.txt"))
  list(networks = networks, audit = audit)
}

compare_aggregation <- function(aoi, pb, out, sensitivity = FALSE) {
  comparisons <- list(); stats <- list()
  for (nm in c("G_DN_vs_Control", "T_DN_vs_Control", "Interaction")) {
    a <- aoi[[nm]]; b <- pb[[nm]]; ii <- match(a$gene, b$gene)
    d <- data.frame(gene = a$gene, AOI_log2FC = a$log2FoldChange, PB_log2FC = b$log2FoldChange[ii],
      AOI_selected = a$Significant, PB_selected = b$Significant[ii])
    ok <- is.finite(d$AOI_log2FC) & is.finite(d$PB_log2FC)
    rho <- stats::cor(d$AOI_log2FC[ok], d$PB_log2FC[ok], method = "spearman")
    d$Selection <- ifelse(d$AOI_selected & d$PB_selected, "Both selected",
      ifelse(d$AOI_selected, "AOI only", ifelse(d$PB_selected, "Pseudobulk only", "Neither selected")))
    comparisons[[nm]] <- d
    stats[[nm]] <- data.frame(Contrast = nm, n = sum(ok), Spearman_rho = rho)
    export <- d
    if (sensitivity) {
      names(export) <- sub("^AOI_", "Primary_", sub("^PB_", "Sensitivity_", names(export)))
      export$Selection <- sub("Pseudobulk only", "Sensitivity only", sub("AOI only", "Primary only", export$Selection, fixed = TRUE), fixed = TRUE)
    }
    write_table(export, file.path(out, paste0(nm, "_effect_comparison.csv")))
  }
  sm <- do.call(rbind, stats); rownames(sm) <- NULL
  write_table(sm, file.path(out, if (sensitivity) "sensitivity_correlations.csv" else "aggregation_correlations.csv"))
  list(tables = comparisons, summary = sm)
}

run_sensitivity <- function(counts, meta, primary, out, workers) {
  pairs <- table(meta$Patient_ID, meta$Compartment)
  complete <- rownames(pairs)[rowSums(pairs > 0) == 2L]
  runs <- c(list(Complete_pairs_only = complete), setNames(
    lapply(levels(meta$Patient_ID), function(x) setdiff(levels(meta$Patient_ID), x)),
    paste0("LOPO_", levels(meta$Patient_ID))))
  logs <- summaries <- list()
  for (nm in names(runs)) {
    keep <- as.character(meta$Patient_ID) %in% runs[[nm]]
    message("Sensitivity: ", nm)
    fit <- tryCatch(fit_dream(counts[, keep, drop = FALSE], meta[keep, , drop = FALSE],
      ~0 + Group + (1 | Patient_ID) + (1 | Patient_Compartment), file.path(out, nm),
      workers, universe = rownames(counts)), error = function(e) e)
    if (inherits(fit, "error")) {
      logs[[nm]] <- data.frame(Run = nm, Status = "FAILED", Detail = conditionMessage(fit))
      next
    }
    x <- classify_interaction(fit$results)
    write_table(x, file.path(out, nm, "Interaction_classified.csv"))
    sm <- compare_aggregation(primary, fit$results, file.path(out, nm, "comparison_to_primary"), sensitivity = TRUE)$summary
    names(sm)[names(sm) == "Spearman_rho"] <- "Spearman_vs_primary"
    sm$Run <- nm; summaries[[nm]] <- sm
    logs[[nm]] <- data.frame(Run = nm, Status = "OK", Detail = "Same model; TMM and voom refitted; unavailable fits retained")
  }
  write_table(do.call(rbind, logs), file.path(out, "sensitivity_run_log.csv"))
  if (length(summaries)) write_table(do.call(rbind, summaries), file.path(out, "sensitivity_effect_summary.csv"))
  if (any(vapply(logs, function(x) x$Status == "FAILED", logical(1))))
    warning("Some sensitivity fits failed; see sensitivity_run_log.csv. No simpler fallback was used.")
}

# --------------------------- Visualization ---------------------------
group_colors <- c(Control_G = "#008000", Control_T = "#8064FF", DN_G = "#FF454B", DN_T = "#DDA51B")
class_colors <- c(G_attenuation = "#0072B2", T_attenuation = "#D55E00",
                  Other_interaction = "#009E73", Not_significant = "#BFBFBF")
disease_colors <- c(Control = "#4477AA", DN = "#CC6677")

save_plot <- function(p, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height,
    device = grDevices::pdf, useDingbats = FALSE, limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height,
    dpi = 300, bg = "white", limitsize = FALSE)
}

volcano <- function(d, contrast) {
  compartment <- contrast %in% c("Control_GT", "DN_GT")
  colors <- switch(contrast,
    Control_GT = c("#8064FF", "#008000"), DN_GT = c("#DDA51B", "#FF454B"),
    G_DN_vs_Control = c("#008000", "#FF454B"), T_DN_vs_Control = c("#8064FF", "#DDA51B"))
  labels <- switch(contrast, Control_GT = c("Higher in Control T", "Higher in Control G"),
    DN_GT = c("Higher in DN T", "Higher in DN G"),
    G_DN_vs_Control = c("Higher in Control G", "Higher in DN G"),
    T_DN_vs_Control = c("Higher in Control T", "Higher in DN T"))
  d$Status <- "Not DEG"
  d$Status[d$Significant & d$log2FoldChange < 0] <- labels[1]
  d$Status[d$Significant & d$log2FoldChange > 0] <- labels[2]
  d <- d[is.finite(d$log2FoldChange) & is.finite(d$padj), , drop = FALSE]
  d$neglogP <- -log10(pmax(d$padj, .Machine$double.xmin))
  labs <- do.call(rbind, lapply(labels, function(x) {
    group <- d[d$Status == x, , drop = FALSE]
    head(group[order(group$padj), , drop = FALSE], 3)
  }))
  xlab <- switch(contrast, Control_GT = "log2 fold change (G / T; Control)",
    DN_GT = "log2 fold change (G / T; DN)", G_DN_vs_Control = "log2 fold change (DN / Control; G)",
    T_DN_vs_Control = "log2 fold change (DN / Control; T)")
  ggplot2::ggplot(d, ggplot2::aes(log2FoldChange, neglogP, color = Status)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.65) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey70") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey70") +
    ggrepel::geom_text_repel(data = labs, ggplot2::aes(label = gene), size = 3, seed = 20261022,
      max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("Not DEG" = "grey80", setNames(colors, labels))) +
    ggplot2::labs(x = xlab, y = "-log10 (BH-adjusted P)", color = NULL) +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}

go_plot <- function(d) {
  if (!nrow(d)) return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
    label = "No eligible displayed GO terms") + ggplot2::theme_void())
  d$Ratio <- vapply(strsplit(d$GeneRatio, "/", fixed = TRUE), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  # Preserve the archived panel order, including ties in GeneRatio.
  d$Term <- factor(d$ID, levels = rev(d$ID))
  lab <- setNames(vapply(d$Description, function(x) paste(strwrap(x, width = 35), collapse = "\n"), character(1)), d$ID)
  ggplot2::ggplot(d, ggplot2::aes(Ratio, Term, size = Count, color = p.adjust)) +
    ggplot2::geom_point(alpha = 0.9) + ggplot2::scale_y_discrete(labels = lab) +
    ggplot2::scale_color_gradient(low = "#EE6666", high = "#2C7FB8", name = "BH-adjusted P") +
    ggplot2::scale_size_continuous(range = c(2, 7)) + ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(x = "GeneRatio", y = NULL) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

long_matrix <- function(x, meta, id_col) {
  d <- data.frame(gene = rep(rownames(x), times = ncol(x)),
    Profile_ID = rep(colnames(x), each = nrow(x)), value = as.vector(x))
  ii <- match(d$Profile_ID, rownames(meta))
  for (nm in c("Group", "Disease", "Patient_ID", "Compartment")) d[[nm]] <- meta[[nm]][ii]
  d$Marker_group <- factor(marker_definitions$Marker_group[match(d$gene, marker_definitions$gene)],
    levels = unique(marker_definitions$Marker_group))
  d$gene <- factor(d$gene, levels = rev(marker_definitions$gene))
  ord <- order(match(as.character(meta$Group), c("Control_G", "Control_T", "DN_G", "DN_T")),
               as.character(meta$Patient_ID), rownames(meta))
  d$Profile_ID <- factor(d$Profile_ID, levels = rownames(meta)[ord])
  d$Group <- factor(d$Group, levels = c("Control_G", "Control_T", "DN_G", "DN_T"))
  d
}

marker_heatmap <- function(x, meta) {
  d <- long_matrix(x, meta, "Profile_ID")
  d$value <- pmax(-3, pmin(3, d$value))
  ggplot2::ggplot(d, ggplot2::aes(Profile_ID, gene, fill = value)) + ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::facet_grid(Marker_group ~ Group, scales = "free", space = "free",
      labeller = ggplot2::labeller(Marker_group = ggplot2::label_wrap_gen(width = 18))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-3, 3), name = "Gene z-score") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      strip.text.y = ggplot2::element_text(angle = 0, size = 9), panel.spacing = grid::unit(0.3, "lines"),
      plot.margin = ggplot2::margin(8, 20, 8, 8)) + ggplot2::labs(x = NULL, y = NULL)
}

make_figures <- function(aoi, pb, meta, profiles, interaction, distance, markers, go, comp, ppi, out) {
  need_packages(c("ggplot2", "ggrepel", "patchwork"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  e <- aoi$logCPM
  v <- apply(e, 1, stats::var)
  genes <- rownames(e)[head(order(v, decreasing = TRUE, na.last = NA), 500)]
  pc <- stats::prcomp(t(e[genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pcvar <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  pd <- data.frame(AOI_ID = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2],
    Group = meta$Group, Disease = meta$Disease)
  write_table(pd, file.path(out, "Figure1A_PCA_coordinates.csv"))
  write_table(data.frame(gene = genes), file.path(out, "Figure1A_PCA_genes.csv"))
  pca <- ggplot2::ggplot(pd, ggplot2::aes(PC1, PC2, color = Group, shape = Disease)) +
    ggplot2::geom_point(size = 2) + ggrepel::geom_text_repel(ggplot2::aes(label = AOI_ID),
      size = 2.4, seed = 20261022, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = group_colors) + ggplot2::scale_shape_manual(values = c(Control = 15, DN = 17)) +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", pcvar[1]), y = sprintf("PC2 (%.1f%%)", pcvar[2]), color = NULL, shape = NULL)
  f1 <- patchwork::wrap_plots(pca, volcano(aoi$results$Control_GT, "Control_GT"),
    volcano(aoi$results$DN_GT, "DN_GT"), nrow = 1) + patchwork::plot_annotation(tag_levels = "A")
  save_plot(f1, file.path(out, "Figure1"), 19, 6)
  selected <- c("CHI3L1", "CLIC5", "ARMH4", "MT1F", "PLG", "SLC39A5")
  x <- interaction[is.finite(interaction$Control_GT_LFC) & is.finite(interaction$DN_GT_LFC), ]
  x <- x[order(x$Class != "Not_significant"), ]
  p2a <- ggplot2::ggplot(x, ggplot2::aes(Control_GT_LFC, DN_GT_LFC, color = Class)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
    ggplot2::geom_hline(yintercept = 0, color = "grey90") + ggplot2::geom_vline(xintercept = 0, color = "grey90") +
    ggplot2::geom_point(size = 1.1, alpha = 0.65) +
    ggrepel::geom_text_repel(data = x[x$gene %in% selected, ], ggplot2::aes(label = gene),
      size = 3, seed = 20261022, show.legend = FALSE, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = class_colors, breaks = names(class_colors),
      labels = c("G attenuation", "T attenuation", "Other interaction", "Not significant")) +
    ggplot2::coord_equal(xlim = c(-5, 5), ylim = c(-5, 5)) +
    ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Control G-T log2 fold change", y = "DN G-T log2 fold change", color = NULL)
  p2b <- ggplot2::ggplot(distance$patient, ggplot2::aes(Disease, RMS_zlogCPM, color = Disease)) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_point(data = distance$summary, ggplot2::aes(Disease, Mean_RMS), inherit.aes = FALSE, shape = 95, size = 7) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Patient_ID), size = 3, seed = 20261022, show.legend = FALSE, max.overlaps = Inf) +
    ggplot2::scale_color_manual(values = disease_colors) + ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none") + ggplot2::labs(x = NULL, y = "Within-patient G-T RMS distance\n(z-scored mean logCPM)")
  xd <- interaction[match(selected, interaction$gene), ]; xd <- xd[!is.na(xd$gene), ]
  xd$gene <- factor(xd$gene, levels = rev(selected))
  p2c <- ggplot2::ggplot(xd, ggplot2::aes(y = gene)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = Control_GT_LFC, xend = DN_GT_LFC, yend = gene), color = "grey65") +
    ggplot2::geom_point(ggplot2::aes(x = Control_GT_LFC, color = "Control"), size = 2.7) +
    ggplot2::geom_point(ggplot2::aes(x = DN_GT_LFC, color = "DN"), size = 2.7) +
    ggplot2::scale_color_manual(values = disease_colors) + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "G-T log2 fold change", y = NULL, color = NULL)
  f2 <- patchwork::wrap_plots(p2a, p2b, p2c, design = "AB\nAC") + patchwork::plot_annotation(tag_levels = "A")
  save_plot(f2, file.path(out, "Figure2"), 16, 9)
  p3a <- volcano(aoi$results$G_DN_vs_Control, "G_DN_vs_Control")
  p3b <- volcano(aoi$results$T_DN_vs_Control, "T_DN_vs_Control")
  if (!is.null(go)) {
    f3 <- patchwork::wrap_plots(p3a, p3b, go_plot(go$displayed$G_Down_in_DN),
      go_plot(go$displayed$T_Down_in_DN), ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
    save_plot(f3, file.path(out, "Figure3"), 17, 15)
    s3 <- patchwork::wrap_plots(go_plot(go$displayed$G_attenuation), go_plot(go$displayed$T_attenuation), ncol = 2) +
      patchwork::plot_annotation(tag_levels = "A")
    save_plot(s3, file.path(out, "FigureS3"), 16, 10)
  } else save_plot(patchwork::wrap_plots(p3a, p3b, ncol = 2) + patchwork::plot_annotation(tag_levels = "A"),
                   file.path(out, "Figure3_AB_only_GO_skipped"), 15, 6)
  save_plot(marker_heatmap(markers$aoi_z, meta), file.path(out, "Figure4"), 18, 11)
  scatterplots <- list()
  cols <- c("Neither selected" = "#BFBFBF", "Pseudobulk only" = "#0072B2", "AOI only" = "#D55E00", "Both selected" = "#009E73")
  for (nm in names(comp$tables)) {
    d <- comp$tables[[nm]]; d <- d[is.finite(d$AOI_log2FC) & is.finite(d$PB_log2FC), ]
    s <- comp$summary[comp$summary$Contrast == nm, ]
    d <- d[order(d$Selection != "Neither selected"), ]
    scatterplots[[nm]] <- ggplot2::ggplot(d, ggplot2::aes(AOI_log2FC, PB_log2FC, color = Selection)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
      ggplot2::geom_point(size = 1, alpha = 0.6) + ggplot2::coord_equal() +
      ggplot2::scale_color_manual(values = cols) + ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "bottom") + ggplot2::labs(title = nm,
        subtitle = sprintf("n = %d; Spearman rho = %.3f", s$n, s$Spearman_rho),
        x = "AOI-dream effect (log2FC)", y = "Pseudobulk-dream effect (log2FC)", color = NULL)
  }
  save_plot(patchwork::wrap_plots(scatterplots, nrow = 1) + patchwork::plot_annotation(tag_levels = "A"),
    file.path(out, "FigureS2"), 18, 6)
  s5a <- marker_heatmap(markers$mean_z, profiles$meta)
  da <- long_matrix(markers$aoi, meta, "AOI_ID")
  dm <- long_matrix(markers$mean, profiles$meta, "Profile_ID")
  da$gene <- factor(da$gene, levels = marker_definitions$gene)
  dm$gene <- factor(dm$gene, levels = marker_definitions$gene)
  s5b <- ggplot2::ggplot(da, ggplot2::aes(Group, value, color = Patient_ID)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.16, height = 0, seed = 20261022), size = 0.9, alpha = 0.6) +
    ggplot2::geom_point(data = dm, shape = 18, size = 1.8) +
    ggplot2::facet_wrap(~gene, scales = "free_y", ncol = 4) +
    ggplot2::scale_x_discrete(labels = c(Control_G = "Control\nG", Control_T = "Control\nT", DN_G = "DN\nG", DN_T = "DN\nT")) +
    ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(size = 6)) +
    ggplot2::labs(x = NULL, y = "TMM-normalized logCPM", color = "Patient")
  save_plot(patchwork::wrap_plots(s5a, s5b, nrow = 1) + patchwork::plot_annotation(tag_levels = "A"),
    file.path(out, "FigureS5"), 21, 10)
  if (!is.null(ppi)) {
    need_packages("igraph")
    graphs <- layouts <- list()
    for (nm in names(ppi$networks)) {
      ed <- ppi$networks[[nm]]
      g <- igraph::graph_from_data_frame(ed[, c("gene_A", "gene_B", "score")], directed = FALSE)
      set.seed(20261022)
      lay <- if (igraph::vcount(g) > 1L) igraph::layout_with_fr(g, niter = 3000) else matrix(c(0, 0), ncol = 2)
      graphs[[nm]] <- g; layouts[[nm]] <- lay
      if (igraph::vcount(g)) write_table(data.frame(gene = igraph::V(g)$name, x = lay[, 1], y = lay[, 2]),
                                        file.path(out, paste0("FigureS4_", nm, "_coordinates.csv")))
    }
    for (ext in c("pdf", "png")) {
      file <- file.path(out, paste0("FigureS4.", ext))
      if (ext == "pdf") grDevices::pdf(file, width = 11, height = 13, useDingbats = FALSE) else
        grDevices::png(file, width = 3300, height = 3900, res = 300)
      graphics::layout(matrix(1:2, ncol = 1), heights = c(5, 1))
      for (nm in names(graphs)) {
        g <- graphs[[nm]]; graphics::par(mar = c(1, 1, 3, 1))
        if (nm == "T_attenuation" && igraph::vcount(g) == 2L && igraph::ecount(g) == 1L) {
          graphics::plot.new()
          graphics::plot.window(xlim = c(-1, 1), ylim = c(-0.5, 0.5))
          graphics::segments(-0.35, 0, 0.35, 0, col = "#89909D", lwd = 2)
          graphics::points(c(-0.35, 0.35), c(0, 0), pch = 21, bg = "#B8E0D2", col = "#688A82", cex = 3)
          graphics::text(c(-0.35, 0.35), c(0.18, 0.18), labels = igraph::V(g)$name, cex = 0.9)
          graphics::title("B  Tubular attenuation")
        } else if (igraph::vcount(g)) {
          graphics::plot(g, layout = layouts[[nm]], vertex.size = if (nm == "G_attenuation") 6 else 18,
            vertex.color = ifelse(igraph::V(g)$name %in% c("SYNPO", "ACTN4"), "#FFD166", "#B8E0D2"),
            vertex.frame.color = "#688A82", vertex.label.color = "black", vertex.label.cex = 0.65,
            edge.width = 1 + 2 * igraph::E(g)$score, edge.color = "#89909D",
            margin = 0.18, main = if (nm == "G_attenuation") "A  Glomerular attenuation" else "B  Tubular attenuation")
        } else { graphics::plot.new(); graphics::title(paste(nm, "- no qualifying edges")) }
      }
      grDevices::dev.off()
    }
  }
}

published_check <- function(aoi, pb, ia, ip, distance, comp, go, ppi, out) {
  rows <- list()
  add <- function(metric, observed, expected, tolerance = 0) {
    rows[[length(rows) + 1L]] <<- data.frame(Metric = metric, Observed = as.numeric(observed),
      Manuscript = expected, Tolerance = tolerance,
      Matches = is.finite(observed) && abs(observed - expected) <= tolerance)
  }
  for (nm in c("G_DN_vs_Control", "T_DN_vs_Control", "Control_GT", "DN_GT")) {
    add(paste0("AOI_", nm, "_selected"), sum(aoi[[nm]]$Significant),
        c(G_DN_vs_Control = 480, T_DN_vs_Control = 141, Control_GT = 484, DN_GT = 58)[nm])
  }
  add("AOI_interaction_testable_genes", sum(is.finite(ia$pvalue)), 18657)
  add("AOI_interactions", sum(ia$Significant), 291)
  add("AOI_G_attenuation", sum(ia$Class == "G_attenuation"), 145)
  add("AOI_T_attenuation", sum(ia$Class == "T_attenuation"), 11)
  add("AOI_other_interaction", sum(ia$Class == "Other_interaction"), 135)
  add("PB_interactions", sum(ip$Significant), 335)
  add("PB_G_attenuation", sum(ip$Class == "G_attenuation"), 166)
  add("PB_T_attenuation", sum(ip$Class == "T_attenuation"), 22)
  for (nm in c("G_DN_vs_Control", "T_DN_vs_Control", "Interaction"))
    add(paste0("Correlation_", nm), comp$summary$Spearman_rho[comp$summary$Contrast == nm],
        c(G_DN_vs_Control = 0.988, T_DN_vs_Control = 0.994, Interaction = 0.979)[nm], 0.0005)
  add("RMS_Control", distance$summary$Mean_RMS[distance$summary$Disease == "Control"], 1.961, 0.0005)
  add("RMS_DN", distance$summary$Mean_RMS[distance$summary$Disease == "DN"], 0.953, 0.0005)
  if (!is.null(go)) for (nm in go$audit$Gene_set) {
    r <- go$audit[go$audit$Gene_set == nm, ]
    add(paste0("GO_", nm, "_effective_foreground"), r$Effective_foreground_n,
        c(G_attenuation = 142, T_attenuation = 10, G_Down_in_DN = 384, T_Down_in_DN = 133)[nm])
    add(paste0("GO_", nm, "_effective_background"), r$Effective_background_n, 16580)
  }
  if (!is.null(ppi)) for (nm in ppi$audit$Gene_set) {
    r <- ppi$audit[ppi$audit$Gene_set == nm, ]
    add(paste0("STRING_", nm, "_connected"), r$Connected_n, c(G_attenuation = 75, T_attenuation = 2)[nm])
    add(paste0("STRING_", nm, "_edges"), r$Edges_n, c(G_attenuation = 113, T_attenuation = 1)[nm])
  }
  d <- do.call(rbind, rows); rownames(d) <- NULL
  write_table(d, file.path(out, "manuscript_comparison.csv"))
  if (any(!d$Matches)) warning("Some results differ from the manuscript. Inspect manuscript_comparison.csv and software/input provenance; results were not forced to match.")
  d
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_options(args)
  if (opt$help) {
    cat(paste(c("Usage: Rscript analysis.R [options]", "",
      "--download       Download the two GEO supplementary inputs (requires public GEO files)",
      "--check-only     Validate counts, AOI identity, cohort and slide assignments; no models",
      "--skip-go        Skip GO enrichment", "--skip-ppi       Skip STRING networks",
      "--skip-plots     Save analysis tables only", "--sensitivity    Also refit complete-pair and leave-one-patient-out AOI models",
      "--string-live    Query STRING instead of using the archived network cache",
      "--workers=N      Parallel model workers; default 1", "--counts=PATH    Override count CSV (or CSV.gz)",
      "--metadata=PATH  Override slide assignment XLSX", "--output=PATH    Empty/new output directory", "--help"), collapse = "\n"), "\n")
    return(invisible(NULL))
  }
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  root <- if (length(script)) dirname(normalizePath(sub("^--file=", "", script[1]))) else normalizePath(getwd())
  data_dir <- file.path(root, "data"); refdir <- file.path(data_dir, "reference")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  if (opt$download) fetch_geo(data_dir)
  count_file <- resolve_input(opt$counts, data_dir, "(^|_)raw_counts\\.csv(\\.gz)?$", "raw_counts.csv")
  meta_file <- resolve_input(opt$metadata, data_dir, "Slide_Assignment\\.xlsx$", "GSE346609_Slide_Assignment.xlsx")
  out <- if (is.null(opt$output)) file.path(root, "results", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())) else opt$output
  if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Output directory must be new or empty: ", out)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  set.seed(20261022)
  on.exit(utils::capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt")), add = TRUE)
  writeLines(c(paste("Started:", Sys.time()), paste("Arguments:", paste(args, collapse = " ")),
    "AOI primary: ~0 + Group + (1|Patient_ID) + (1|Patient_Compartment)",
    "Pseudobulk: ~0 + Group + (1|Patient_ID)",
    "TMM; prior.count=0.5; rescaleWeightsAfter=FALSE; REML; Satterthwaite; eBayes",
    "FDR<0.05 per contrast; abs(log2FC)>1 additionally for non-interaction DEGs",
    "No AOI exclusion, predefined-gene exclusion or additional LOQ filtering",
    "QC: only area/nuclei flags can be reconstructed from the deposited slide workbook",
    "Selected-gene RMS is descriptive; no confirmatory P value",
    "Slide is not fitted as a separate effect; disease and slide are confounded"), file.path(out, "analysis_settings.txt"))
  inputs <- read_inputs(count_file, meta_file)
  counts <- inputs$counts; meta <- inputs$meta
  write_table(meta, file.path(out, "metadata_AOI.csv"))
  write_table(inputs$filter, file.path(out, "gene_filter_audit.csv"))
  write_table(data.frame(File = basename(c(count_file, meta_file)), MD5 = unname(tools::md5sum(c(count_file, meta_file)))),
              file.path(out, "input_checksums.csv"))
  slide <- as.data.frame(table(Slide = meta$SlideName, Disease = meta$Disease))
  write_table(slide, file.path(out, "slide_by_disease.csv"))
  mm <- model.matrix(~Disease + factor(SlideName), meta)
  rank_info <- data.frame(Design = "Disease + Slide", Columns = ncol(mm), Rank = qr(mm)$rank,
                         Separately_estimable = qr(mm)$rank == ncol(mm))
  write_table(rank_info, file.path(out, "slide_confounding_audit.csv"))
  write_table(as.data.frame(table(Group = meta$Group)), file.path(out, "AOI_group_counts.csv"))
  writeLines(c(paste("Genes retained:", nrow(counts)), paste("AOIs:", ncol(counts)),
    paste("Patients:", nlevels(meta$Patient_ID)), paste("Nuclei <200:", sum(meta$Low_nuclei_flag)),
    "Input AOI IDs are matched by name; Normal is standardized to Control.",
    "No missing glomerular profile is imputed."), file.path(out, "input_validation.txt"))
  message("Validated ", nrow(counts), " genes, ", ncol(counts), " AOIs and ", nlevels(meta$Patient_ID), " patients.")
  if (opt$check_only) { message("Input validation complete: ", normalizePath(out)); return(invisible(inputs)) }
  pkgs <- c("edgeR", "limma", "variancePartition", "BiocParallel")
  if (!opt$skip_go) pkgs <- c(pkgs, "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
  if (!opt$skip_plots) pkgs <- c(pkgs, "ggplot2", "ggrepel", "patchwork")
  if (!opt$skip_ppi && !opt$skip_plots) pkgs <- c(pkgs, "igraph")
  if (!opt$skip_ppi && opt$string_live) pkgs <- c(pkgs, "rbioapi")
  need_packages(pkgs)
  if (utils::packageVersion("variancePartition") < "1.32.5") stop("variancePartition >=1.32.5 is required.")
  record_versions(pkgs, refdir, out)
  message("Fitting the primary AOI model.")
  aoi <- fit_dream(counts, meta, ~0 + Group + (1 | Patient_ID) + (1 | Patient_Compartment), file.path(out, "AOI"), opt$workers)
  ia <- classify_interaction(aoi$results)
  write_table(ia, file.path(out, "AOI", "Interaction_classified.csv"))
  for (cls in c("G_attenuation", "T_attenuation", "Other_interaction"))
    write_table(ia[ia$Class == cls, , drop = FALSE], file.path(out, "AOI", paste0(cls, ".csv")))
  profiles <- aggregate_profiles(counts, aoi$logCPM, meta)
  write_matrix(profiles$counts, file.path(out, "Pseudobulk", "summed_raw_counts.csv"))
  write_matrix(profiles$means, file.path(out, "Display", "patient_compartment_mean_logCPM.csv"))
  write_table(profiles$meta, file.path(out, "Display", "patient_compartment_metadata.csv"))
  message("Fitting the 15-profile patient-compartment pseudobulk model.")
  pb <- fit_dream(profiles$counts, profiles$meta, ~0 + Group + (1 | Patient_ID), file.path(out, "Pseudobulk"), opt$workers)
  ip <- classify_interaction(pb$results)
  write_table(ip, file.path(out, "Pseudobulk", "Interaction_classified.csv"))
  comp <- compare_aggregation(aoi$results, pb$results, file.path(out, "Aggregation"))
  distance <- selected_distance(profiles$means, profiles$meta, aoi$results, file.path(out, "Distance"))
  markers <- prepare_markers(aoi$logCPM, meta, profiles$means, profiles$meta, file.path(out, "Markers"))
  for (direction in c("Down_in_DN", "Up_in_DN")) {
    r <- aoi$results
    shared <- intersect(r$G_DN_vs_Control$gene[r$G_DN_vs_Control$Direction == direction],
                        r$T_DN_vs_Control$gene[r$T_DN_vs_Control$Direction == direction])
    write_table(data.frame(gene = shared), file.path(out, "AOI", paste0("Shared_", direction, ".csv")))
  }
  go <- if (opt$skip_go) NULL else run_go(aoi$results, ia, file.path(out, "GO"), refdir)
  ppi <- if (opt$skip_ppi) NULL else run_ppi(ia, file.path(out, "STRING"), refdir, opt$string_live)
  check <- published_check(aoi$results, pb$results, ia, ip, distance, comp, go, ppi, out)
  if (!opt$skip_plots) make_figures(aoi, pb, meta, profiles, ia, distance, markers, go, comp, ppi, file.path(out, "Figures"))
  if (opt$sensitivity) run_sensitivity(counts, meta, aoi$results, file.path(out, "Sensitivity"), opt$workers)
  writeLines(c(paste("Completed:", Sys.time()), paste("GO performed:", !opt$skip_go),
    paste("STRING performed:", !opt$skip_ppi), paste("Plots generated:", !opt$skip_plots),
    paste("Optional sensitivity requested:", opt$sensitivity), paste("All checked manuscript quantities match:", all(check$Matches))),
    file.path(out, "RUN_COMPLETE.txt"))
  message("Analysis complete: ", normalizePath(out))
  invisible(list(output = out, manuscript_comparison = check))
}

# Setting this environment variable allows isolated function tests without
# running the full analysis. Normal Rscript/source use runs main().
if (Sys.getenv("DN_RELEASE_SOURCE_ONLY") != "1") main()
