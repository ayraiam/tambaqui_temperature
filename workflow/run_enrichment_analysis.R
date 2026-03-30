#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Dr.eg.db)
  library(enrichplot)
})

message2 <- function(...) cat(..., "\n", sep = "")

ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

save_table <- function(x, path) {
  if (is.null(x)) return(invisible(NULL))
  write.table(as.data.frame(x), file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

save_dotplot <- function(obj, path, title_txt, show_n = 20) {
  if (is.null(obj)) return(invisible(NULL))
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(invisible(NULL))
  p <- dotplot(obj, showCategory = min(show_n, nrow(df))) + ggtitle(title_txt)
  ggsave(path, plot = p, width = 8, height = 6, dpi = 300)
}

dedup_ranked <- function(df, gene_col, score_col) {
  df |>
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[score_col]])) |>
    group_by(.data[[gene_col]]) |>
    slice_max(order_by = abs(.data[[score_col]]), n = 1, with_ties = FALSE) |>
    ungroup()
}

split_core_enrichment <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  trimws(unlist(strsplit(x, "/", fixed = TRUE)))
}

count_gsea_core_terms <- function(gsea_df, alpha = 0.05, target_geneid_col = "target_Gene_ID") {
  gsea_sig <- gsea_df |>
    dplyr::filter(!is.na(p.adjust), p.adjust <= alpha) |>
    dplyr::filter(!is.na(core_enrichment), core_enrichment != "")
  
  if (nrow(gsea_sig) == 0) {
    out <- tibble::tibble(tmp_id = character(), gsea_go_core_term_count = integer())
    colnames(out)[1] <- target_geneid_col
    return(out)
  }
  
  core_long <- gsea_sig |>
    dplyr::rowwise() |>
    dplyr::mutate(tmp_id = list(split_core_enrichment(core_enrichment))) |>
    tidyr::unnest(cols = c(tmp_id)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(tmp_id), tmp_id != "")
  
  out <- core_long |>
    dplyr::count(tmp_id, name = "gsea_go_core_term_count") |>
    dplyr::arrange(dplyr::desc(gsea_go_core_term_count), tmp_id)
  
  colnames(out)[1] <- target_geneid_col
  out
}

make_candidate_gene_table <- function(
    annotated_tsv,
    normalized_counts_tsv,
    metadata_tsv,
    gsea_go_tsv,
    outdir,
    sample_col = "sample",
    group_col = "Condition",
    alpha = 0.05,
    source_geneid_col = "Geneid",
    target_geneid_col = "target_Gene_ID"
) {
  message2(">>> Building candidate gene table")
  
  annotated_df <- read.delim(annotated_tsv, check.names = FALSE)
  norm_mat <- read.delim(normalized_counts_tsv, check.names = FALSE, row.names = 1)
  meta_df <- read.delim(metadata_tsv, check.names = FALSE)
  gsea_df <- read.delim(gsea_go_tsv, check.names = FALSE)
  
  norm_mat <- as.matrix(norm_mat)
  
  required_annotated <- c(source_geneid_col, target_geneid_col, "log2FoldChange", "padj", "baseMean")
  missing_annotated <- setdiff(required_annotated, colnames(annotated_df))
  if (length(missing_annotated) > 0) {
    stop("Annotated table is missing required columns: ",
         paste(missing_annotated, collapse = ", "))
  }
  
  if (!(sample_col %in% colnames(meta_df))) {
    stop("Metadata file is missing sample column: ", sample_col)
  }
  
  if (!(group_col %in% colnames(meta_df))) {
    stop("Metadata file is missing group column: ", group_col)
  }
  
  if (!("core_enrichment" %in% colnames(gsea_df))) {
    stop("GSEA GO table is missing core_enrichment column")
  }
  
  if (!("GeneName" %in% colnames(annotated_df))) {
    annotated_df$GeneName <- NA_character_
  }
  
  meta_df[[sample_col]] <- as.character(meta_df[[sample_col]])
  meta_df[[group_col]] <- as.character(meta_df[[group_col]])
  
  common_samples <- meta_df[[sample_col]][meta_df[[sample_col]] %in% colnames(norm_mat)]
  common_samples <- unique(common_samples)
  
  if (length(common_samples) == 0) {
    stop("No overlapping samples between metadata and normalized counts matrix")
  }
  
  meta_df <- meta_df |>
    dplyr::filter(.data[[sample_col]] %in% common_samples)
  
  norm_mat <- norm_mat[, common_samples, drop = FALSE]
  
  groups <- unique(meta_df[[group_col]])
  if (length(groups) != 2) {
    stop("Candidate table step expects exactly 2 groups in metadata_used.tsv for column '",
         group_col, "'. Found: ", paste(groups, collapse = ", "))
  }
  
  group1 <- groups[1]
  group2 <- groups[2]
  
  group1_samples <- meta_df |>
    dplyr::filter(.data[[group_col]] == group1) |>
    dplyr::pull(.data[[sample_col]])
  
  group2_samples <- meta_df |>
    dplyr::filter(.data[[group_col]] == group2) |>
    dplyr::pull(.data[[sample_col]])
  
  mean_norm_all_df <- tibble::tibble(
    join_source_geneid = rownames(norm_mat),
    mean_norm_all = rowMeans(norm_mat, na.rm = TRUE)
  )
  
  mean_norm_g1_df <- tibble::tibble(
    join_source_geneid = rownames(norm_mat),
    value = rowMeans(norm_mat[, group1_samples, drop = FALSE], na.rm = TRUE)
  )
  colnames(mean_norm_g1_df)[2] <- paste0("mean_norm_", group1)
  
  mean_norm_g2_df <- tibble::tibble(
    join_source_geneid = rownames(norm_mat),
    value = rowMeans(norm_mat[, group2_samples, drop = FALSE], na.rm = TRUE)
  )
  colnames(mean_norm_g2_df)[2] <- paste0("mean_norm_", group2)
  
  gsea_counts_df <- count_gsea_core_terms(
    gsea_df,
    alpha = alpha,
    target_geneid_col = target_geneid_col
  )
  
  annotated_df[[source_geneid_col]] <- as.character(annotated_df[[source_geneid_col]])
  annotated_df[[target_geneid_col]] <- as.character(annotated_df[[target_geneid_col]])
  
  candidate_df <- annotated_df |>
    dplyr::mutate(
      abs_log2FoldChange = abs(log2FoldChange)
    ) |>
    dplyr::left_join(
      mean_norm_all_df,
      by = stats::setNames("join_source_geneid", source_geneid_col)
    ) |>
    dplyr::left_join(
      mean_norm_g1_df,
      by = stats::setNames("join_source_geneid", source_geneid_col)
    ) |>
    dplyr::left_join(
      mean_norm_g2_df,
      by = stats::setNames("join_source_geneid", source_geneid_col)
    ) |>
    dplyr::left_join(
      gsea_counts_df,
      by = target_geneid_col
    ) |>
    dplyr::mutate(
      gsea_go_core_term_count = ifelse(is.na(gsea_go_core_term_count), 0L, gsea_go_core_term_count),
      candidate_score = abs_log2FoldChange *
        log10(mean_norm_all + 1) *
        (1 + log2(gsea_go_core_term_count + 1))
    )
  
  candidate_df <- candidate_df |>
    dplyr::rename(Geneid = !!source_geneid_col) |>
    dplyr::select(
      Geneid,
      GeneName,
      log2FoldChange,
      abs_log2FoldChange,
      padj,
      baseMean,
      mean_norm_all,
      dplyr::starts_with("mean_norm_"),
      gsea_go_core_term_count,
      candidate_score
    ) |>
    dplyr::arrange(
      dplyr::desc(gsea_go_core_term_count),
      dplyr::desc(abs_log2FoldChange),
      padj,
      dplyr::desc(baseMean)
    )
  
  out_tsv <- file.path(outdir, "candidate_gene_table.tsv")
  write.table(
    candidate_df,
    file = out_tsv,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  summary_df <- tibble::tibble(
    metric = c(
      "n_rows_in_candidate_table",
      "n_genes_with_gsea_core_support",
      "group1",
      "group2",
      "n_group1_samples",
      "n_group2_samples"
    ),
    value = c(
      nrow(candidate_df),
      sum(candidate_df$gsea_go_core_term_count > 0, na.rm = TRUE),
      group1,
      group2,
      length(group1_samples),
      length(group2_samples)
    )
  )
  
  write.table(
    summary_df,
    file = file.path(outdir, "candidate_gene_table_summary.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  message2(">>> Candidate gene table written to: ", out_tsv)
}

get_args <- function() {
  option_list <- list(
    make_option("--task", type = "character", dest = "task", default = "analysis",
                help = "Task to run: analysis or candidates [default: %default]"),
    make_option("--annotated-tsv", type = "character", dest = "annotated_tsv",
                help = "Annotated DE table with Danio rerio mappings"),
    make_option("--normalized-counts-tsv", type = "character", dest = "normalized_counts_tsv",
                help = "DESeq2 normalized_counts.tsv"),
    make_option("--metadata-tsv", type = "character", dest = "metadata_tsv",
                help = "DESeq2 metadata_used.tsv"),
    make_option("--gsea-go-tsv", type = "character", dest = "gsea_go_tsv",
                help = "GSEA GO BP TSV"),
    make_option("--sample-col", type = "character", dest = "sample_col", default = "sample",
                help = "Sample column in metadata [default: %default]"),
    make_option("--group-col", type = "character", dest = "group_col", default = "Condition",
                help = "Grouping column in metadata [default: %default]"),
    make_option("--outdir", type = "character", dest = "outdir",
                help = "Output directory"),
    make_option("--alpha", type = "double", dest = "alpha", default = 0.05,
                help = "Adjusted p-value cutoff [default: %default]"),
    make_option("--target-geneid-col", type = "character", dest = "target_geneid_col", default = "target_Gene_ID",
                help = "Column containing zebrafish Gene IDs [default: %default]"),
    make_option("--source-geneid-col", type = "character", dest = "source_geneid_col", default = "Geneid",
                help = "Source organism gene ID column [default: %default]"),
    make_option("--logfc-col", type = "character", dest = "logfc_col", default = "log2FoldChange",
                help = "Column containing ranking statistic / log2FC [default: %default]"),
    make_option("--padj-col", type = "character", dest = "padj_col", default = "padj",
                help = "Column containing adjusted p-values [default: %default]")
  )
  
  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)
  
  if (args$task == "analysis") {
    required <- c("annotated_tsv", "outdir")
  } else if (args$task == "candidates") {
    required <- c("annotated_tsv", "normalized_counts_tsv", "metadata_tsv", "gsea_go_tsv", "outdir")
  } else {
    print_help(parser)
    stop("Unsupported --task value: ", args$task)
  }
  
  missing_required <- required[vapply(required, function(x) is.null(args[[x]]), logical(1))]
  if (length(missing_required) > 0) {
    print_help(parser)
    stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
  }
  
  args
}

run_enrichment_main <- function(args) {
  df <- read.delim(args$annotated_tsv, check.names = FALSE)
  
  stopifnot(args$target_geneid_col %in% colnames(df))
  stopifnot(args$logfc_col %in% colnames(df))
  stopifnot(args$padj_col %in% colnames(df))
  
  df[[args$target_geneid_col]] <- as.character(df[[args$target_geneid_col]])
  df[[args$logfc_col]] <- as.numeric(df[[args$logfc_col]])
  df[[args$padj_col]] <- as.numeric(df[[args$padj_col]])
  
  mapped_all <- df |>
    filter(!is.na(.data[[args$target_geneid_col]]), !is.na(.data[[args$logfc_col]]))
  
  universe_ids <- unique(mapped_all[[args$target_geneid_col]])
  
  sig <- mapped_all |>
    filter(!is.na(.data[[args$padj_col]]), .data[[args$padj_col]] <= args$alpha)
  
  up_ids <- unique(sig |>
                     filter(.data[[args$logfc_col]] > 0) |>
                     pull(.data[[args$target_geneid_col]]))
  
  down_ids <- unique(sig |>
                       filter(.data[[args$logfc_col]] < 0) |>
                       pull(.data[[args$target_geneid_col]]))
  
  ranked_df <- dedup_ranked(mapped_all, args$target_geneid_col, args$logfc_col)
  ranked <- ranked_df[[args$logfc_col]]
  names(ranked) <- ranked_df[[args$target_geneid_col]]
  ranked <- sort(ranked, decreasing = TRUE)
  
  write.table(
    data.frame(target_Gene_ID = names(ranked), ranking_score = ranked),
    file = file.path(args$outdir, "gsea_ranked_gene_list.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  ora_go_up <- if (length(up_ids) > 0) enrichGO(
    gene = up_ids,
    universe = universe_ids,
    OrgDb = org.Dr.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = args$alpha,
    qvalueCutoff = args$alpha,
    readable = TRUE
  ) else NULL
  
  ora_go_down <- if (length(down_ids) > 0) enrichGO(
    gene = down_ids,
    universe = universe_ids,
    OrgDb = org.Dr.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = args$alpha,
    qvalueCutoff = args$alpha,
    readable = TRUE
  ) else NULL
  
  ora_kegg_up <- if (length(up_ids) > 0) enrichKEGG(
    gene = up_ids,
    universe = universe_ids,
    organism = "dre",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = args$alpha,
    qvalueCutoff = args$alpha
  ) else NULL
  
  ora_kegg_down <- if (length(down_ids) > 0) enrichKEGG(
    gene = down_ids,
    universe = universe_ids,
    organism = "dre",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = args$alpha,
    qvalueCutoff = args$alpha
  ) else NULL
  
  save_table(ora_go_up, file.path(args$outdir, "ora_go_bp_up.tsv"))
  save_table(ora_go_down, file.path(args$outdir, "ora_go_bp_down.tsv"))
  save_table(ora_kegg_up, file.path(args$outdir, "ora_kegg_up.tsv"))
  save_table(ora_kegg_down, file.path(args$outdir, "ora_kegg_down.tsv"))
  
  save_dotplot(ora_go_up, file.path(args$outdir, "ora_go_bp_up_dotplot.png"), "ORA GO BP - up")
  save_dotplot(ora_go_down, file.path(args$outdir, "ora_go_bp_down_dotplot.png"), "ORA GO BP - down")
  save_dotplot(ora_kegg_up, file.path(args$outdir, "ora_kegg_up_dotplot.png"), "ORA KEGG - up")
  save_dotplot(ora_kegg_down, file.path(args$outdir, "ora_kegg_down_dotplot.png"), "ORA KEGG - down")
  
  gsea_go <- if (length(ranked) > 1) gseGO(
    geneList = ranked,
    OrgDb = org.Dr.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = args$alpha,
    verbose = FALSE
  ) else NULL
  
  gsea_kegg <- if (length(ranked) > 1) gseKEGG(
    geneList = ranked,
    organism = "dre",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = args$alpha,
    verbose = FALSE
  ) else NULL
  
  save_table(gsea_go, file.path(args$outdir, "gsea_go_bp.tsv"))
  save_table(gsea_kegg, file.path(args$outdir, "gsea_kegg.tsv"))
  
  save_dotplot(gsea_go, file.path(args$outdir, "gsea_go_bp_dotplot.png"), "GSEA GO BP")
  save_dotplot(gsea_kegg, file.path(args$outdir, "gsea_kegg_dotplot.png"), "GSEA KEGG")
  
  summary_df <- data.frame(
    metric = c(
      "mapped_rows",
      "universe_ids",
      "sig_rows",
      "up_ids",
      "down_ids"
    ),
    value = c(
      nrow(mapped_all),
      length(universe_ids),
      nrow(sig),
      length(up_ids),
      length(down_ids)
    )
  )
  write.table(summary_df,
              file = file.path(args$outdir, "enrichment_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message2(">>> Enrichment outputs written to: ", args$outdir)
}

main <- function() {
  args <- get_args()
  ensure_dir(args$outdir)
  
  if (args$task == "analysis") {
    run_enrichment_main(args)
  } else if (args$task == "candidates") {
    make_candidate_gene_table(
      annotated_tsv = args$annotated_tsv,
      normalized_counts_tsv = args$normalized_counts_tsv,
      metadata_tsv = args$metadata_tsv,
      gsea_go_tsv = args$gsea_go_tsv,
      outdir = args$outdir,
      sample_col = args$sample_col,
      group_col = args$group_col,
      alpha = args$alpha,
      source_geneid_col = args$source_geneid_col,
      target_geneid_col = args$target_geneid_col
    )
  }
}

main()
