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

get_args <- function() {
  option_list <- list(
    make_option("--annotated-tsv", type = "character", dest = "annotated_tsv",
                help = "Annotated DE table with Danio rerio mappings"),
    make_option("--outdir", type = "character", dest = "outdir",
                help = "Output directory"),
    make_option("--alpha", type = "double", dest = "alpha", default = 0.05,
                help = "Adjusted p-value cutoff [default: %default]"),
    make_option("--target-geneid-col", type = "character", dest = "target_geneid_col", default = "target_Gene_ID",
                help = "Column containing zebrafish Gene IDs [default: %default]"),
    make_option("--logfc-col", type = "character", dest = "logfc_col", default = "log2FoldChange",
                help = "Column containing ranking statistic / log2FC [default: %default]"),
    make_option("--padj-col", type = "character", dest = "padj_col", default = "padj",
                help = "Column containing adjusted p-values [default: %default]")
  )
  
  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)
  
  required <- c("annotated_tsv", "outdir")
  missing_required <- required[vapply(required, function(x) is.null(args[[x]]), logical(1))]
  
  if (length(missing_required) > 0) {
    print_help(parser)
    stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
  }
  
  args
}

main <- function() {
  args <- get_args()
  ensure_dir(args$outdir)
  
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
  
  # ORA
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
  
  # GSEA
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

main()
