#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(edgeR)
  library(limma)
  library(variancePartition)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript workflow/plot_variance_partition.R <featureCounts_tsv> <metadata_tsv> <outdir>")
}

counts_tsv <- args[[1]]
metadata_tsv <- args[[2]]
outdir <- args[[3]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Settings
# -----------------------------
excluded_samples <- c("RFA-64", "RFA-70", "RFA-74")

extract_sample_id <- function(x) {
  x2 <- ifelse(grepl("/", x, fixed = TRUE), basename(dirname(x)), x)
  sub("^(RFA-[0-9]+).*", "\\1", x2)
}

exclude_flagged_samples <- function(x) {
  x[!x %in% excluded_samples]
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

# -----------------------------
# Load counts
# -----------------------------
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

# Remove excluded samples
keep_samples <- exclude_flagged_samples(colnames(count_matrix))
count_matrix <- count_matrix[, keep_samples, drop = FALSE]

# Order samples
sample_levels <- ordered_rfa_levels(colnames(count_matrix))
count_matrix <- count_matrix[, sample_levels, drop = FALSE]

# -----------------------------
# Load metadata
# -----------------------------
meta <- read_tsv(metadata_tsv, show_col_types = FALSE) %>%
  mutate(sample = extract_sample_id(sample)) %>%
  filter(!sample %in% excluded_samples)

# Match order
meta <- meta %>%
  filter(sample %in% colnames(count_matrix)) %>%
  mutate(sample = factor(sample, levels = colnames(count_matrix))) %>%
  arrange(sample)

# -----------------------------
# Build Age and Temperature
# -----------------------------
# Expected metadata columns:
# - condition (C0, T1, T2)
# - age (1, 2 OR "1 month", "2 months")
# - temperature (pre, room, 37C)

# Clean Age
meta$Age <- as.character(meta$age)
meta$Age[meta$Age == "1"] <- "1 month"
meta$Age[meta$Age == "2"] <- "2 months"
meta$Age <- factor(meta$Age, levels = c("1 month", "2 months"))

# Clean Temperature
meta$Temperature <- as.character(meta$temperature)
meta$Temperature <- factor(meta$Temperature)

# -----------------------------
# edgeR + voom
# -----------------------------
dge <- DGEList(counts = count_matrix)

# Filter low expression genes
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ 1, data = meta)

v <- voom(dge, design, plot = FALSE)

# -----------------------------
# Variance Partition
# -----------------------------
form <- ~ Age + Temperature

varPart <- fitExtractVarPartModel(v$E, form, meta)

# -----------------------------
# Save gene-level results
# -----------------------------
varpart_df <- as.data.frame(varPart)
varpart_df$gene <- rownames(varpart_df)

write_tsv(
  varpart_df,
  file.path(outdir, "variance_partition_gene_level.tsv")
)

# -----------------------------
# Summary (median variance)
# -----------------------------
summary_df <- varpart_df %>%
  select(-gene) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(cols = everything(), names_to = "factor", values_to = "median_variance")

write_tsv(
  summary_df,
  file.path(outdir, "variance_partition_summary.tsv")
)

# -----------------------------
# Prepare plot
# -----------------------------
plot_df <- varpart_df %>%
  select(-gene) %>%
  pivot_longer(cols = everything(), names_to = "factor", values_to = "variance")

plot_df$factor <- recode(
  plot_df$factor,
  Age = "Age",
  Temperature = "Temperature",
  Residuals = "Residual"
)

plot_df$factor <- factor(
  plot_df$factor,
  levels = c("Age", "Temperature", "Residual")
)

# -----------------------------
# Plot (violin + boxplot)
# -----------------------------
p <- ggplot(plot_df, aes(x = factor, y = variance, fill = factor)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c(
    "Age" = "#8da0cb",
    "Temperature" = "#fc8d62",
    "Residual" = "#bdbdbd"
  )) +
  labs(
    title = "",
    x = NULL,
    y = "Fraction of variance explained\n"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(outdir, "variance_partition_violin.pdf"),
  p,
  width = 5,
  height = 5
)

ggsave(
  file.path(outdir, "variance_partition_violin.png"),
  p,
  width = 5,
  height = 5,
  dpi = 300
)

message(">>> Variance partition analysis complete.")
message(">>> Outputs written to: ", outdir)
