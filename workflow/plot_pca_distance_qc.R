#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript workflow/plot_pca_distance_qc.R <featureCounts_tsv> <metadata_tsv> <outdir>")
}

counts_tsv <- args[[1]]
metadata_tsv <- args[[2]]
outdir <- args[[3]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

excluded_samples <- c("RFA-64", "RFA-70", "RFA-74")
base_font_size <- 16

exclude_flagged_samples <- function(x) {
  x[!x %in% excluded_samples]
}

extract_sample_id <- function(x) {
  x2 <- ifelse(grepl("/", x, fixed = TRUE), basename(dirname(x)), x)
  sub("^(RFA-[0-9]+).*", "\\1", x2)
}

ordered_rfa_levels <- function(x) {
  tibble(sample = unique(x)) %>%
    mutate(
      sample = extract_sample_id(sample),
      sample_num = as.numeric(sub("^RFA-([0-9]+).*", "\\1", sample))
    ) %>%
    arrange(sample_num, sample) %>%
    pull(sample)
}

get_condition_palette <- function(conditions) {
  base_pal <- c(
    "T1" = "darkmagenta",
    "C0" = "olivedrab",
    "T2" = "goldenrod2"
  )
  
  conds <- unique(as.character(conditions))
  missing <- setdiff(conds, names(base_pal))
  
  if (length(missing) > 0) {
    extra_cols <- setNames(RColorBrewer::brewer.pal(length(missing), "Set2"), missing)
    base_pal <- c(base_pal, extra_cols)
  }
  
  base_pal[conds]
}

load_featurecounts_matrix <- function(counts_tsv) {
  counts <- read.delim(
    counts_tsv,
    comment.char = "#",
    check.names = FALSE
  )
  
  count_matrix <- counts[, 7:ncol(counts), drop = FALSE]
  rownames(count_matrix) <- make.unique(as.character(counts$Geneid))
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  
  sample_names <- colnames(count_matrix)
  sample_names <- extract_sample_id(sample_names)
  colnames(count_matrix) <- sample_names
  
  if (anyDuplicated(colnames(count_matrix))) {
    stop(
      "Duplicate sample names detected after parsing featureCounts columns: ",
      paste(unique(colnames(count_matrix)[duplicated(colnames(count_matrix))]), collapse = ", ")
    )
  }
  
  keep_samples <- exclude_flagged_samples(colnames(count_matrix))
  count_matrix <- count_matrix[, keep_samples, drop = FALSE]
  
  count_matrix
}

load_metadata <- function(metadata_tsv, sample_names) {
  meta <- read_tsv(metadata_tsv, show_col_types = FALSE) %>%
    mutate(sample = extract_sample_id(sample)) %>%
    filter(!sample %in% excluded_samples)
  
  required_cols <- c("sample", "condition")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (anyDuplicated(meta$sample)) {
    stop(
      "Duplicate sample IDs found in metadata: ",
      paste(unique(meta$sample[duplicated(meta$sample)]), collapse = ", ")
    )
  }
  
  missing_in_meta <- setdiff(sample_names, meta$sample)
  if (length(missing_in_meta) > 0) {
    stop(
      "These samples are in featureCounts but missing from metadata: ",
      paste(missing_in_meta, collapse = ", ")
    )
  }
  
  meta <- meta %>%
    filter(sample %in% sample_names) %>%
    mutate(sample = factor(sample, levels = sample_names)) %>%
    arrange(sample)
  
  if ("age" %in% colnames(meta)) {
    meta$age <- as.character(meta$age)
    meta$age[meta$age == "1"] <- "1 month"
    meta$age[meta$age == "2"] <- "2 months"
    meta$age <- factor(meta$age, levels = c("1 month", "2 months"))
  }
  
  as.data.frame(meta)
}

make_vst <- function(count_matrix, metadata_df) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = metadata_df,
    design = ~ condition
  )
  vst(dds, blind = TRUE)
}

plot_pca_from_vst <- function(vsd, metadata_df, out_prefix, outdir, title_text) {
  pca <- prcomp(t(assay(vsd)))
  percent_var <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100
  
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$sample <- rownames(scores)
  scores <- left_join(scores, metadata_df, by = "sample")
  
  write_tsv(scores, file.path(outdir, paste0(out_prefix, "_scores.tsv")))
  
  cond_pal <- get_condition_palette(scores$condition)
  
  p <- ggplot(scores, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5.5) +
    scale_color_manual(values = cond_pal) +
    labs(
      title = title_text,
      x = paste0("PC1 (", round(percent_var[1], 1), "%)\n"),
      y = paste0("\nPC2 (", round(percent_var[2], 1), "%)"),
      color = "Condition"
    ) +
    theme_classic(base_size = base_font_size) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(out_prefix, ".pdf")),
    plot = p,
    width = 6,
    height = 5
  )
  
  ggsave(
    filename = file.path(outdir, paste0(out_prefix, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

plot_pca_global_and_t1t2 <- function(count_matrix, metadata_df, outdir) {
  vsd_all <- make_vst(count_matrix, metadata_df)
  plot_pca_from_vst(
    vsd = vsd_all,
    metadata_df = metadata_df,
    out_prefix = "pca_global_all_samples",
    outdir = outdir,
    title_text = ""
  )
  
  t1t2_samples <- metadata_df$condition %in% c("T1", "T2")
  if (sum(t1t2_samples) < 2) {
    stop("Not enough T1/T2 samples found in metadata to build T1 vs T2 PCA")
  }
  
  metadata_t1t2 <- metadata_df[t1t2_samples, , drop = FALSE]
  count_matrix_t1t2 <- count_matrix[, metadata_t1t2$sample, drop = FALSE]
  vsd_t1t2 <- make_vst(count_matrix_t1t2, metadata_t1t2)
  
  plot_pca_from_vst(
    vsd = vsd_t1t2,
    metadata_df = metadata_t1t2,
    out_prefix = "pca_t1_vs_t2",
    outdir = outdir,
    title_text = ""
  )
}

plot_sample_distance_heatmap <- function(count_matrix, metadata_df, outdir) {
  vsd <- make_vst(count_matrix, metadata_df)
  
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  rownames(sample_dist_matrix) <- metadata_df$sample
  colnames(sample_dist_matrix) <- metadata_df$sample
  
  annotation_df <- metadata_df
  rownames(annotation_df) <- annotation_df$sample
  annotation_df$sample <- NULL
  
  colnames(annotation_df) <- dplyr::recode(
    colnames(annotation_df),
    condition = "Condition",
    age = "Age",
    temperature = "Temperature"
  )
  
  ann_colors <- list()
  
  if ("Condition" %in% colnames(annotation_df)) {
    ann_colors$Condition <- get_condition_palette(annotation_df$Condition)
  }
  
  if ("Age" %in% colnames(annotation_df)) {
    age_levels <- levels(annotation_df$Age)
    ann_colors$Age <- c(
      "1 month" = "#F0027F",
      "2 months" = "#FFFF99"
    )[age_levels]
  }
  
  if ("Temperature" %in% colnames(annotation_df)) {
    temp_levels <- unique(as.character(annotation_df$Temperature))
    temp_palette <- c(
      "pre" = "#7fc97f",
      "room" = "#beaed4",
      "37C" = "#fdc086"
    )
    ann_colors$Temperature <- temp_palette[temp_levels]
  }
  
  heat_cols <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  write_tsv(
    data.frame(sample = rownames(sample_dist_matrix), sample_dist_matrix, check.names = FALSE),
    file.path(outdir, "sample_distance_matrix_all_samples.tsv")
  )
  
  pdf(file.path(outdir, "sample_distance_heatmap_all_samples.pdf"), width = 8, height = 7)
  pheatmap(
    sample_dist_matrix,
    color = heat_cols,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = ann_colors,
    fontsize = 14,
    border_color = NA,
    main = ""
  )
  dev.off()
  
  png(file.path(outdir, "sample_distance_heatmap_all_samples.png"), width = 2600, height = 2500, res = 300)
  pheatmap(
    sample_dist_matrix,
    color = heat_cols,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = ann_colors,
    fontsize = 14,
    border_color = NA,
    main = ""
  )
  dev.off()
}

count_matrix <- load_featurecounts_matrix(counts_tsv)
sample_levels <- ordered_rfa_levels(colnames(count_matrix))
count_matrix <- count_matrix[, sample_levels, drop = FALSE]

metadata_df <- load_metadata(metadata_tsv, colnames(count_matrix))

plot_pca_global_and_t1t2(count_matrix, metadata_df, outdir)
plot_sample_distance_heatmap(count_matrix, metadata_df, outdir)

message(">>> Wrote PCA and sample-distance heatmap outputs to: ", outdir)
