#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(DESeq2)
})

# ==========================================================
# Helpers
# ==========================================================

message2 <- function(...) cat(..., "\n", sep = "")

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

clean_sample_names <- function(x) {
  x |>
    basename() |>
    str_replace("\\.bam$", "") |>
    str_replace("\\.sortedByCoord\\.out$", "") |>
    str_replace("\\.Aligned$", "") |>
    str_replace("\\.txt$", "") |>
    str_replace("\\.tsv$", "")
}

parse_csv_arg <- function(x) {
  if (is.null(x) || is.na(x) || x == "") return(character(0))
  trimws(unlist(strsplit(x, ",")))
}

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(filename, plot = plot_obj, width = width, height = height, dpi = 300)
}

# ==========================================================
# Argument parsing
# ==========================================================

get_args <- function() {
  option_list <- list(
    make_option("--counts-tsv", type = "character", dest = "counts_tsv",
                help = "featureCounts TSV file"),
    make_option("--metadata-tsv", type = "character", dest = "metadata_tsv",
                help = "Metadata TSV file"),
    make_option("--outdir", type = "character", dest = "outdir",
                help = "Output directory"),
    
    make_option("--sample-col", type = "character", dest = "sample_col", default = "sample",
                help = "Sample column in metadata [default: %default]"),
    
    make_option("--design", type = "character", dest = "design", default = "~ Condition",
                help = "DESeq2 design formula [default: %default]"),
    
    make_option("--subset-column", type = "character", dest = "subset_column", default = NULL,
                help = "Optional metadata column used to subset samples"),
    make_option("--subset-values", type = "character", dest = "subset_values", default = NULL,
                help = "Comma-separated values to keep from subset-column"),
    
    make_option("--reference-variable", type = "character", dest = "reference_variable", default = NULL,
                help = "Factor column to relevel before DESeq2"),
    make_option("--reference-level", type = "character", dest = "reference_level", default = NULL,
                help = "Reference level for reference-variable"),
    
    make_option("--contrast-variable", type = "character", dest = "contrast_variable", default = NULL,
                help = "Variable for DESeq2 contrast, e.g. Condition"),
    make_option("--contrast-numerator", type = "character", dest = "contrast_numerator", default = NULL,
                help = "Numerator level, e.g. T2"),
    make_option("--contrast-denominator", type = "character", dest = "contrast_denominator", default = NULL,
                help = "Denominator/reference level, e.g. T1"),
    
    make_option("--alpha", type = "double", dest = "alpha", default = 0.05,
                help = "Adjusted p-value cutoff [default: %default]"),
    
    make_option("--min-count", type = "integer", dest = "min_count", default = 10,
                help = "Minimum count threshold for gene filtering [default: %default]"),
    make_option("--min-samples", type = "integer", dest = "min_samples", default = 3,
                help = "Minimum number of samples meeting min-count [default: %default]"),
    
    make_option("--vst-blind", action = "store_true", dest = "vst_blind", default = FALSE,
                help = "Run vst(blind=TRUE). Default is FALSE"),
    
    make_option("--annotation-tsv", type = "character", dest = "annotation_tsv", default = NULL,
                help = "Optional gene annotation TSV"),
    make_option("--annotation-id-col", type = "character", dest = "annotation_id_col", default = "Geneid",
                help = "Gene ID column in annotation TSV [default: %default]"),
    make_option("--annotation-name-col", type = "character", dest = "annotation_name_col", default = "GeneName",
                help = "Gene symbol/name column in annotation TSV [default: %default]")
  )
  
  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)
  print(args)
  
  required <- c("counts_tsv", "metadata_tsv", "outdir")
  missing_required <- required[vapply(required, function(x) is.null(args[[x]]), logical(1))]
  
  if (length(missing_required) > 0) {
    print_help(parser)
    stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
  }
  
  args
}

# ==========================================================
# Input reading
# ==========================================================

read_featurecounts <- function(counts_tsv) {
  message2(">>> Reading featureCounts file: ", counts_tsv)
  
  fc <- read.delim(counts_tsv, comment.char = "#", check.names = FALSE)
  
  required_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  missing_cols <- setdiff(required_cols, colnames(fc))
  if (length(missing_cols) > 0) {
    stop("featureCounts file is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  sample_cols <- setdiff(colnames(fc), required_cols)
  if (length(sample_cols) == 0) {
    stop("No sample count columns found in featureCounts file.")
  }
  
  counts_mat <- fc |>
    dplyr::select(all_of(sample_cols)) |>
    as.matrix()
  
  mode(counts_mat) <- "integer"
  rownames(counts_mat) <- fc$Geneid
  colnames(counts_mat) <- clean_sample_names(colnames(counts_mat))
  
  gene_info <- fc |>
    dplyr::select(all_of(required_cols))
  
  list(
    counts = counts_mat,
    gene_info = gene_info
  )
}

read_metadata_table <- function(metadata_tsv, sample_col) {
  message2(">>> Reading metadata file: ", metadata_tsv)
  
  meta <- read.delim(metadata_tsv, check.names = FALSE)
  
  if (!(sample_col %in% colnames(meta))) {
    stop("Metadata file does not contain sample column: ", sample_col)
  }
  
  meta[[sample_col]] <- clean_sample_names(meta[[sample_col]])
  meta <- meta |>
    distinct(.data[[sample_col]], .keep_all = TRUE)
  
  meta
}

subset_metadata <- function(meta, subset_column = NULL, subset_values = NULL) {
  if (is.null(subset_column) || is.null(subset_values)) {
    return(meta)
  }
  
  keep_values <- parse_csv_arg(subset_values)
  
  if (!(subset_column %in% colnames(meta))) {
    stop("subset-column not found in metadata: ", subset_column)
  }
  
  meta2 <- meta |>
    filter(.data[[subset_column]] %in% keep_values)
  
  if (nrow(meta2) == 0) {
    stop("No samples remaining after subsetting on ", subset_column,
         " with values: ", paste(keep_values, collapse = ", "))
  }
  
  meta2
}

align_counts_and_metadata <- function(counts_mat, meta, sample_col) {
  common_samples <- intersect(colnames(counts_mat), meta[[sample_col]])
  
  if (length(common_samples) == 0) {
    stop("No overlapping sample names between counts matrix and metadata.")
  }
  
  counts_mat2 <- counts_mat[, common_samples, drop = FALSE]
  meta2 <- meta[match(common_samples, meta[[sample_col]]), , drop = FALSE]
  
  rownames(meta2) <- meta2[[sample_col]]
  
  if (!identical(colnames(counts_mat2), rownames(meta2))) {
    stop("Counts and metadata could not be aligned correctly.")
  }
  
  list(
    counts = counts_mat2,
    meta = meta2
  )
}

apply_reference_level <- function(meta, reference_variable = NULL, reference_level = NULL) {
  if (is.null(reference_variable) || is.null(reference_level)) {
    return(meta)
  }
  
  if (!(reference_variable %in% colnames(meta))) {
    stop("reference-variable not found in metadata: ", reference_variable)
  }
  
  meta[[reference_variable]] <- factor(meta[[reference_variable]])
  
  if (!(reference_level %in% levels(meta[[reference_variable]]))) {
    stop("reference-level '", reference_level, "' not found in metadata column '",
         reference_variable, "'.")
  }
  
  meta[[reference_variable]] <- relevel(meta[[reference_variable]], ref = reference_level)
  meta
}

coerce_design_variables <- function(meta, design_formula) {
  design_vars <- all.vars(as.formula(design_formula))
  
  for (v in design_vars) {
    if (!(v %in% colnames(meta))) {
      stop("Design variable not found in metadata: ", v)
    }
    
    if (is.character(meta[[v]]) || is.logical(meta[[v]])) {
      meta[[v]] <- factor(meta[[v]])
    }
  }
  
  meta
}

# ==========================================================
# DESeq2 setup and filtering
# ==========================================================

build_dds <- function(counts_mat, meta, design_formula) {
  message2(">>> Building DESeqDataSet with design: ", design_formula)
  
  DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData = meta,
    design = as.formula(design_formula)
  )
}

filter_low_count_genes <- function(dds, min_count = 10, min_samples = 3, outdir = NULL) {
  message2(">>> Filtering low-count genes")
  n_before <- nrow(dds)
  
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds_filtered <- dds[keep, ]
  
  n_after <- nrow(dds_filtered)
  
  stats_df <- data.frame(
    step = c("before_filter", "after_filter"),
    n_genes = c(n_before, n_after)
  )
  
  if (!is.null(outdir)) {
    write.table(
      stats_df,
      file = file.path(outdir, "filtering_summary.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  
  message2("    Genes before filtering: ", n_before)
  message2("    Genes after filtering : ", n_after)
  
  dds_filtered
}

# ==========================================================
# Transformed data for downstream visualization
# ==========================================================

run_vst_transform <- function(dds, blind = FALSE) {
  message2(">>> Running VST transformation (blind = ", blind, ")")
  vst(dds, blind = blind)
}

save_transformed_and_normalized_counts <- function(dds, vsd, outdir) {
  norm_counts <- counts(dds, normalized = TRUE)
  vst_mat <- assay(vsd)
  
  write.table(
    norm_counts,
    file = file.path(outdir, "normalized_counts.tsv"),
    sep = "\t", quote = FALSE, col.names = NA
  )
  
  write.table(
    vst_mat,
    file = file.path(outdir, "vst_matrix.tsv"),
    sep = "\t", quote = FALSE, col.names = NA
  )
}

# ==========================================================
# Differential expression
# ==========================================================

run_deseq_model <- function(dds) {
  message2(">>> Running DESeq()")
  DESeq(dds)
}

extract_deseq_results <- function(dds,
                                  contrast_variable,
                                  contrast_numerator,
                                  contrast_denominator,
                                  alpha = 0.05) {
  if (is.null(contrast_variable) ||
      is.null(contrast_numerator) ||
      is.null(contrast_denominator)) {
    stop("Contrast arguments are required: --contrast-variable, --contrast-numerator, --contrast-denominator")
  }
  
  contrast_vec <- c(contrast_variable, contrast_numerator, contrast_denominator)
  
  message2(">>> Extracting results for contrast: ",
           contrast_variable, " | ",
           contrast_numerator, " vs ",
           contrast_denominator)
  
  res <- results(
    dds,
    contrast = contrast_vec,
    alpha = alpha,
    independentFiltering = TRUE,
    pAdjustMethod = "BH"
  )
  
  shrink_ok <- requireNamespace("apeglm", quietly = TRUE)
  
  if (shrink_ok) {
    res_shrunk <- lfcShrink(
      dds,
      contrast = contrast_vec,
      res = res,
      type = "apeglm"
    )
  } else {
    message2(">>> apeglm not available; returning unshrunk log2FC")
    res_shrunk <- res
  }
  
  as.data.frame(res_shrunk) |>
    rownames_to_column("Geneid")
}

annotate_results <- function(res_df,
                             annotation_tsv = NULL,
                             annotation_id_col = "Geneid",
                             annotation_name_col = "GeneName") {
  if (is.null(annotation_tsv)) {
    return(res_df)
  }
  
  ann <- read.delim(annotation_tsv, check.names = FALSE)
  
  if (!(annotation_id_col %in% colnames(ann))) {
    stop("annotation-id-col not found in annotation TSV: ", annotation_id_col)
  }
  if (!(annotation_name_col %in% colnames(ann))) {
    stop("annotation-name-col not found in annotation TSV: ", annotation_name_col)
  }
  
  ann2 <- ann |>
    dplyr::select(all_of(annotation_id_col), all_of(annotation_name_col)) |>
    distinct()
  
  colnames(ann2) <- c("Geneid", "GeneName")
  
  left_join(res_df, ann2, by = "Geneid")
}

save_result_tables <- function(res_df, alpha, outdir, prefix) {
  res_df <- res_df |>
    mutate(
      significant = ifelse(!is.na(padj) & padj <= alpha, "yes", "no"),
      direction = case_when(
        is.na(log2FoldChange) ~ NA_character_,
        log2FoldChange > 0 ~ "up",
        log2FoldChange < 0 ~ "down",
        TRUE ~ "flat"
      )
    ) |>
    arrange(padj)
  
  sig_df <- res_df |>
    filter(!is.na(padj), padj <= alpha)
  
  write.table(
    res_df,
    file = file.path(outdir, paste0(prefix, "_all_genes.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  write.table(
    sig_df,
    file = file.path(outdir, paste0(prefix, "_significant_genes.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  summary_df <- data.frame(
    metric = c("total_genes_tested", "significant_genes", "upregulated", "downregulated"),
    value = c(
      sum(!is.na(res_df$padj)),
      nrow(sig_df),
      sum(sig_df$log2FoldChange > 0, na.rm = TRUE),
      sum(sig_df$log2FoldChange < 0, na.rm = TRUE)
    )
  )
  
  write.table(
    summary_df,
    file = file.path(outdir, paste0(prefix, "_summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  list(all = res_df, sig = sig_df, summary = summary_df)
}

plot_ma_result <- function(res_df, outdir, prefix, alpha = 0.05) {
  plot_df <- res_df |>
    filter(!is.na(baseMean), !is.na(log2FoldChange)) |>
    mutate(sig = !is.na(padj) & padj <= alpha)
  
  p <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(
      title = "MA-style plot",
      x = "log10(baseMean + 1)",
      y = "log2 fold change"
    )
  
  save_plot(p, file.path(outdir, paste0(prefix, "_MA_plot.png")), width = 7, height = 5)
}

plot_volcano_result <- function(res_df, outdir, prefix, alpha = 0.05) {
  plot_df <- res_df |>
    mutate(
      neglog10padj = -log10(padj),
      status = case_when(
        is.na(padj) ~ "NS",
        padj <= alpha & log2FoldChange > 0 ~ "Up",
        padj <= alpha & log2FoldChange < 0 ~ "Down",
        padj <= alpha & log2FoldChange == 0 ~ "Sig",
        TRUE ~ "NS"
      )
    )
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10padj, color = status)) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(
      title = "Volcano plot",
      x = "log2 fold change",
      y = "-log10 adjusted p-value"
    )
  
  save_plot(p, file.path(outdir, paste0(prefix, "_volcano.png")), width = 7, height = 5)
}

save_session_info <- function(outdir) {
  writeLines(capture.output(sessionInfo()),
             con = file.path(outdir, "sessionInfo.txt"))
}

# ==========================================================
# Main
# ==========================================================

main <- function() {
  args <- get_args()
  ensure_dir(args$outdir)
  
  fc <- read_featurecounts(args$counts_tsv)
  meta <- read_metadata_table(args$metadata_tsv, args$sample_col)
  meta <- subset_metadata(meta, args$subset_column, args$subset_values)
  meta <- apply_reference_level(meta, args$reference_variable, args$reference_level)
  meta <- coerce_design_variables(meta, args$design)
  
  aligned <- align_counts_and_metadata(fc$counts, meta, args$sample_col)
  
  counts_mat <- aligned$counts
  meta_aligned <- aligned$meta
  
  write.table(
    meta_aligned,
    file = file.path(args$outdir, "metadata_used.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  dds <- build_dds(counts_mat, meta_aligned, args$design)
  
  dds <- filter_low_count_genes(
    dds,
    min_count = args$min_count,
    min_samples = args$min_samples,
    outdir = args$outdir
  )
  
  vsd <- run_vst_transform(dds, blind = args$vst_blind)
  save_transformed_and_normalized_counts(dds, vsd, args$outdir)
  
  dds <- run_deseq_model(dds)
  
  res_df <- extract_deseq_results(
    dds = dds,
    contrast_variable = args$contrast_variable,
    contrast_numerator = args$contrast_numerator,
    contrast_denominator = args$contrast_denominator,
    alpha = args$alpha
  )
  
  res_df <- annotate_results(
    res_df,
    annotation_tsv = args$annotation_tsv,
    annotation_id_col = args$annotation_id_col,
    annotation_name_col = args$annotation_name_col
  )
  
  prefix <- paste(
    args$contrast_variable,
    args$contrast_numerator,
    "vs",
    args$contrast_denominator,
    sep = "_"
  )
  
  save_result_tables(
    res_df = res_df,
    alpha = args$alpha,
    outdir = args$outdir,
    prefix = prefix
  )
  
  plot_ma_result(
    res_df = res_df,
    outdir = args$outdir,
    prefix = prefix,
    alpha = args$alpha
  )
  
  plot_volcano_result(
    res_df = res_df,
    outdir = args$outdir,
    prefix = prefix,
    alpha = args$alpha
  )
  
  save_session_info(args$outdir)
  
  message2(">>> Done.")
  message2(">>> Results written to: ", args$outdir)
}

main()
