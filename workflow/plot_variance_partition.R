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

missing_in_meta <- setdiff(colnames(count_matrix), meta$sample)
if (length(missing_in_meta) > 0) {
  stop(
    "These samples are in featureCounts but missing from metadata: ",
    paste(missing_in_meta, collapse = ", ")
  )
}

meta <- meta %>%
  filter(sample %in% colnames(count_matrix)) %>%
  mutate(
    sample = as.character(sample),
    Condition = factor(condition, levels = c("C0", "T1", "T2"))
  )

meta <- as.data.frame(meta)
rownames(meta) <- meta$sample
meta <- meta[colnames(count_matrix), , drop = FALSE]

if (!identical(rownames(meta), colnames(count_matrix))) {
  stop("Metadata rownames do not match count matrix column names.")
}

# -----------------------------
# edgeR + voom
# -----------------------------
dge <- DGEList(counts = count_matrix)

keep <- filterByExpr(dge, group = meta$Condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ Condition, data = meta)
v <- voom(dge, design, plot = FALSE)

# -----------------------------
# Variance Partition
# -----------------------------
form <- ~ Condition
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
  Condition = "Condition",
  Residuals = "Residual"
)

plot_df$factor <- factor(
  plot_df$factor,
  levels = c("Condition", "Residual")
)

# -----------------------------
# Plot (violin + boxplot)
# -----------------------------
p <- ggplot(plot_df, aes(x = factor, y = variance, fill = factor)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c(
    "Condition" = "grey", #"#8da0cb"
    "Residual" = "grey",#"#bdbdbd"
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
  width = 2.5,
  height = 5,
  dpi = 300
)

message(">>> Variance partition analysis complete.")
message(">>> Outputs written to: ", outdir)
