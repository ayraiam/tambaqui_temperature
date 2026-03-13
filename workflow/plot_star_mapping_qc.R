#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript workflow/plot_star_mapping_qc.R <summary_tsv> <outdir> <featureCounts_tsv>")
}

summary_tsv <- args[[1]]
outdir <- args[[2]]
counts_tsv <- args[[3]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

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
    mutate(
      sample = sub("^(RFA-[0-9]+).*", "\\1", sample)
    )
  
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
    scale_fill_brewer(palette = "Paired") +
    labs(
      title = "",
      x = NULL,
      y = "\nReads (%)",
      fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      ),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, "star_mapping_qc_stacked_percent.pdf"),
    plot = p_pct,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = file.path(outdir, "star_mapping_qc_stacked_percent.png"),
    plot = p_pct,
    width = 10,
    height = 7,
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
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8),
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
  sample_names <- basename(sample_names)
  sample_names <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", sample_names)
  sample_names <- sub("\\.bam$", "", sample_names)
  sample_names <- sub("^(RFA-[0-9]+).*", "\\1", sample_names)
  colnames(count_matrix) <- sample_names
  
  library_sizes <- colSums(count_matrix)
  detected_genes <- colSums(count_matrix > 0)
  
  qc_table <- data.frame(
    sample = colnames(count_matrix),
    library_size = as.numeric(library_sizes),
    detected_genes = as.numeric(detected_genes),
    stringsAsFactors = FALSE
  ) %>%
    distinct() %>%
    mutate(
      sample_num = as.numeric(sub("^RFA-([0-9]+).*", "\\1", sample))
    ) %>%
    arrange(sample_num, sample) %>%
    select(-sample_num)
  
  sample_levels <- qc_table$sample
  
  write_tsv(qc_table, file.path(outdir, "featurecounts_sample_qc.tsv"))
  
  p_lib <- ggplot(
    qc_table,
    aes(x = factor(sample, levels = sample_levels), y = library_size)
  ) +
    geom_col(width = 0.8) +
    labs(
      title = "",
      x = NULL,
      y = "\nAssigned reads"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      )
    )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_library_size_per_sample.pdf"),
    plot = p_lib,
    width = 10,
    height = 5
  )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_library_size_per_sample.png"),
    plot = p_lib,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  p_det <- ggplot(
    qc_table,
    aes(x = factor(sample, levels = sample_levels), y = detected_genes)
  ) +
    geom_col(width = 0.8) +
    labs(
      title = "",
      x = NULL,
      y = "\nDetected genes"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      )
    )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_detected_genes_per_sample.pdf"),
    plot = p_det,
    width = 10,
    height = 5
  )
  
  ggsave(
    filename = file.path(outdir, "featurecounts_detected_genes_per_sample.png"),
    plot = p_det,
    width = 10,
    height = 5,
    dpi = 300
  )
}

plot_mapping_qc(summary_tsv, outdir)
plot_featurecounts_sample_qc(counts_tsv, outdir)

message(">>> Wrote mapping QC and featureCounts sample QC outputs to: ", outdir)
