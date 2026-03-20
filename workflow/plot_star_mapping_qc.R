#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript workflow/plot_star_mapping_qc.R <summary_tsv> <outdir> <featureCounts_tsv>")
}

summary_tsv <- args[[1]]
outdir <- args[[2]]
counts_tsv <- args[[3]]

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

sample_dot_palette <- c(
  "#c31f22",
  "#e5ebf5",
  "#000000",
  "#7b7b7b",
  "#8ecdcb",
  "#1e6289",
  "#efe444",
  "#e69f05",
  "#7c1718",
  "#4a4a4a",
  "#bfc7d5",
  "#4fb6a8",
  "#2f8fb3",
  "#f4a261",
  "#ffd166"
)

mapping_qc_palette <- c(
  "Uniquely mapped" = "#a6cee3",
  "Mapped to multiple loci" = "#1f78b4",
  "Mapped to too many loci" = "#b2df8a",
  "Unmapped: too many mismatches" = "#fb9a99",
  "Unmapped: too short" = "#e31a1c",
  "Unmapped: other" = "#fdbf6f",
  "Chimeric" = "#ff7f00"
)

ordered_rfa_levels <- function(x) {
  tibble(sample = unique(x)) %>%
    mutate(
      sample = sub("^(RFA-[0-9]+).*", "\\1", sample),
      sample_num = as.numeric(sub("^RFA-([0-9]+).*", "\\1", sample))
    ) %>%
    arrange(sample_num, sample) %>%
    pull(sample)
}

plot_mapping_qc <- function(summary_tsv, outdir) {
  df <- read_tsv(summary_tsv, show_col_types = FALSE) %>%
    mutate(sample = extract_sample_id(sample)) %>%
    filter(!sample %in% excluded_samples)
  
  sample_levels <- ordered_rfa_levels(df$sample)
  
  pct_long <- df %>%
    select(
      sample,
      uniquely_mapped_reads_pct,
      mapped_multiple_loci_pct,
      mapped_too_many_loci_pct,
      unmapped_too_many_mismatches_pct,
      unmapped_too_short_pct,
      unmapped_other_pct,
      chimeric_reads_pct
    ) %>%
    pivot_longer(
      cols = -sample,
      names_to = "category",
      values_to = "percent"
    ) %>%
    mutate(
      category = recode(
        category,
        uniquely_mapped_reads_pct = "Uniquely mapped",
        mapped_multiple_loci_pct = "Mapped to multiple loci",
        mapped_too_many_loci_pct = "Mapped to too many loci",
        unmapped_too_many_mismatches_pct = "Unmapped: too many mismatches",
        unmapped_too_short_pct = "Unmapped: too short",
        unmapped_other_pct = "Unmapped: other",
        chimeric_reads_pct = "Chimeric"
      ),
      sample = factor(sample, levels = sample_levels)
    )
  
  cnt_long <- df %>%
    select(
      sample,
      uniquely_mapped_reads_n,
      mapped_multiple_loci_n,
      mapped_too_many_loci_n,
      unmapped_too_many_mismatches_n,
      unmapped_too_short_n,
      unmapped_other_n,
      chimeric_reads_n
    ) %>%
    pivot_longer(
      cols = -sample,
      names_to = "category",
      values_to = "count"
    ) %>%
    mutate(
      category = recode(
        category,
        uniquely_mapped_reads_n = "Uniquely mapped",
        mapped_multiple_loci_n = "Mapped to multiple loci",
        mapped_too_many_loci_n = "Mapped to too many loci",
        unmapped_too_many_mismatches_n = "Unmapped: too many mismatches",
        unmapped_too_short_n = "Unmapped: too short",
        unmapped_other_n = "Unmapped: other",
        chimeric_reads_n = "Chimeric"
      ),
      sample = factor(sample, levels = sample_levels)
    )
  
  write_tsv(pct_long, file.path(outdir, "star_mapping_qc_percent_long.tsv"))
  write_tsv(cnt_long, file.path(outdir, "star_mapping_qc_counts_long.tsv"))
  
  p_pct <- ggplot(pct_long, aes(x = sample, y = percent, fill = category)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = mapping_qc_palette) +
    labs(
      title = "",
      x = NULL,
      y = "\nReads (%)",
      fill = NULL
    ) +
    theme_classic(base_size = base_font_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 13
      ),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, "star_mapping_qc_stacked_percent.pdf"),
    plot = p_pct,
    width = 10,
    height = 12
  )
  
  ggsave(
    filename = file.path(outdir, "star_mapping_qc_stacked_percent.png"),
    plot = p_pct,
    width = 10,
    height = 12,
    dpi = 300
  )
  
  p_cnt <- ggplot(cnt_long, aes(x = sample, y = count, fill = category)) +
    geom_col(width = 0.8) +
    coord_flip() +
    labs(
      title = "STAR mapping QC composition by sample",
      x = NULL,
      y = "Reads (count)",
      fill = NULL
    ) +
    theme_bw(base_size = base_font_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 13),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, "star_mapping_qc_stacked_counts.pdf"),
    plot = p_cnt,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = file.path(outdir, "star_mapping_qc_stacked_counts.png"),
    plot = p_cnt,
    width = 10,
    height = 7,
    dpi = 300
  )
}

plot_featurecounts_sample_qc <- function(counts_tsv, outdir) {
  counts <- read.delim(
    counts_tsv,
    comment.char = "#",
    check.names = FALSE
  )
  
  count_matrix <- counts[, 7:ncol(counts), drop = FALSE]
  rownames(count_matrix) <- counts$Geneid
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  
  sample_names <- colnames(count_matrix)
  sample_names <- basename(dirname(sample_names))
  sample_names <- sub("^(RFA-[0-9]+).*", "\\1", sample_names)
  colnames(count_matrix) <- sample_names
  
  if (anyDuplicated(colnames(count_matrix))) {
    stop(
      "Duplicate sample names detected after parsing featureCounts columns: ",
      paste(unique(colnames(count_matrix)[duplicated(colnames(count_matrix))]), collapse = ", ")
    )
  }
  
  keep_samples <- exclude_flagged_samples(colnames(count_matrix))
  count_matrix <- count_matrix[, keep_samples, drop = FALSE]
  
  library_sizes <- colSums(count_matrix)
  detected_genes <- colSums(count_matrix > 0)
  
  qc_table <- data.frame(
    sample = colnames(count_matrix),
    library_size = as.numeric(library_sizes),
    detected_genes = as.numeric(detected_genes),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      sample_num = as.numeric(sub("^RFA-([0-9]+).*", "\\1", sample))
    ) %>%
    arrange(sample_num, sample) %>%
    select(-sample_num)
  
  sample_levels <- qc_table$sample
  
  write_tsv(qc_table, file.path(outdir, "featurecounts_sample_qc.tsv"))
  
  set.seed(123)
  dot_colors <- sample(sample_dot_palette, size = nrow(qc_table), replace = FALSE)
  
  qc_table_lib <- qc_table %>%
    mutate(
      sample = factor(sample, levels = sample_levels),
      dot_color = dot_colors,
      group = "Assigned reads"
    )
  
  p_lib <- ggplot(
    qc_table_lib,
    aes(x = group, y = library_size)
  ) +
    geom_violin(
      fill = "grey90",
      color = "transparent",
      width = 0.9,
      alpha = 0.75,
      trim = FALSE
    ) +
    geom_beeswarm(
      aes(color = sample),
      size = 3.2,
      cex = 3.5,
      priority = "density"
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = NA,
      color = "black",
      alpha = .75
    ) +
    coord_flip() +
    scale_color_manual(values = setNames(qc_table_lib$dot_color, qc_table_lib$sample)) +
    labs(
      title = "",
      x = NULL,
      y = "Assigned reads\n"
    ) +
    theme_classic(base_size = base_font_size) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_library_size_violin_box_beeswarm.pdf"),
    plot = p_lib,
    width = 5,
    height = 4.5
  )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_library_size_violin_box_beeswarm.png"),
    plot = p_lib,
    width = 5,
    height = 4.5,
    dpi = 300
  )
  
  qc_table_det <- qc_table %>%
    mutate(
      sample = factor(sample, levels = sample_levels),
      dot_color = dot_colors,
      group = "Detected genes"
    )
  
  p_det <- ggplot(
    qc_table_det,
    aes(x = group, y = detected_genes)
  ) +
    geom_violin(
      fill = "grey90",
      color = "transparent",
      width = 0.9,
      alpha = 0.75,
      trim = FALSE
    ) +
    geom_beeswarm(
      aes(color = sample),
      size = 3.2,
      cex = 3.5,
      priority = "density"
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = NA,
      color = "black",
      alpha = .75
    ) +
    coord_flip() +
    scale_color_manual(values = setNames(qc_table_det$dot_color, qc_table_det$sample)) +
    labs(
      title = "",
      x = NULL,
      y = "Detected genes\n"
    ) +
    theme_classic(base_size = base_font_size) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_detected_genes_violin_box_beeswarm.pdf"),
    plot = p_det,
    width = 5,
    height = 4.5
  )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_detected_genes_violin_box_beeswarm.png"),
    plot = p_det,
    width = 5,
    height = 4.5,
    dpi = 300
  )
}

plot_mapping_qc(summary_tsv, outdir)
plot_featurecounts_sample_qc(counts_tsv, outdir)

message(">>> Wrote mapping QC and featureCounts sample QC outputs to: ", outdir)
