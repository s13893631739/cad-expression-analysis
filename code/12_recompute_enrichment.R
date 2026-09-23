# Recompute GO and KEGG enrichment from the project differential-expression run.
# Outputs are written to the current project output directory.

script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

experiment_root <- Sys.getenv(
  "CAD_EXPERIMENT_ROOT",
  unset = file.path(output_root, "20_full_experiment")
)
diff_root <- file.path(experiment_root, "03_diff")
out_root <- file.path(output_root, "12_enrichment_recomputed")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

required <- file.path(diff_root, c("Batch.diff.txt", "Batch.all.txt"))
if (any(!file.exists(required))) {
  stop(
    "Differential-expression outputs are missing. Run `Rscript code/run_project.R` first: ",
    paste(required[!file.exists(required)], collapse = ", ")
  )
}

diff_table <- read.delim(required[1], stringsAsFactors = FALSE, check.names = FALSE)
all_table <- read.delim(required[2], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("id", "P.Value") %in% names(diff_table)) || !"id" %in% names(all_table)) {
  stop("Differential-expression tables must contain an `id` column.")
}

deg_symbols <- unique(diff_table$id[nzchar(diff_table$id)])
universe_symbols <- unique(all_table$id[nzchar(all_table$id)])

gene_map <- bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
gene_map <- unique(gene_map[, c("SYMBOL", "ENTREZID")])
deg_map <- gene_map[gene_map$SYMBOL %in% deg_symbols, , drop = FALSE]
if (nrow(deg_map) == 0 || nrow(gene_map) == 0) {
  stop("No genes could be mapped from SYMBOL to ENTREZID.")
}

genes <- unique(deg_map$ENTREZID)
universe <- unique(gene_map$ENTREZID)

params <- c(
  paste("experiment_root", experiment_root),
  paste("de_gene_symbols", length(deg_symbols)),
  paste("universe_gene_symbols", length(universe_symbols)),
  paste("mapped_de_entrez", length(genes)),
  paste("mapped_universe_entrez", length(universe)),
  "organism hsa",
  "kegg_source KEGG REST API (use_internal_data FALSE)",
  "annotation_orgdb org.Hs.eg.db",
  paste("orgdb_version", as.character(packageVersion("org.Hs.eg.db"))),
  paste("clusterprofiler_version", as.character(packageVersion("clusterProfiler"))),
  "minGSSize 10",
  "maxGSSize 500",
  "pAdjustMethod BH",
  "all_results_pvalueCutoff 1",
  "all_results_qvalueCutoff 1",
  "significant_results_p_adjusted < 0.05",
  "qvalue is recorded but is not used as an additional filter; significance is p.adjust < 0.05"
)
writeLines(params, file.path(out_root, "PARAMETERS.txt"))
write.table(gene_map, file.path(out_root, "symbol_entrez_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(deg_map, file.path(out_root, "differential_symbol_entrez_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

as_result_table <- function(x) {
  if (is.null(x)) return(data.frame())
  as.data.frame(x, stringsAsFactors = FALSE)
}

write_result_table <- function(df, path, columns) {
  if (nrow(df) == 0) {
    # Keep a readable header even when an enrichment service returns no rows
    # or is unavailable; a blank file is ambiguous to downstream readers.
    df <- as.data.frame(setNames(replicate(length(columns), character(), simplify = FALSE), columns),
                        stringsAsFactors = FALSE)
  }
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

go <- enrichGO(
  gene = genes,
  universe = universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)
go_df <- as_result_table(go)
if (nrow(go_df) > 0) {
  go_sig <- go_df[is.finite(go_df$p.adjust) & go_df$p.adjust < 0.05, , drop = FALSE]
} else {
  go_sig <- go_df
}
go_columns <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
                "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
write_result_table(go_df, file.path(out_root, "go_all.tsv"), go_columns)
write_result_table(go_sig, file.path(out_root, "go_significant.tsv"), go_columns)

kegg_status <- "completed"
kegg_error <- character()
kegg <- tryCatch(
  enrichKEGG(
    gene = genes,
    universe = universe,
    organism = "hsa",
    keyType = "ncbi-geneid",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    use_internal_data = FALSE
  ),
  error = function(e) {
    kegg_status <<- "failed"
    kegg_error <<- conditionMessage(e)
    NULL
  }
)
kegg_df <- as_result_table(kegg)
if (nrow(kegg_df) > 0) {
  kegg_sig <- kegg_df[is.finite(kegg_df$p.adjust) & kegg_df$p.adjust < 0.05, , drop = FALSE]
} else {
  kegg_sig <- kegg_df
}
kegg_columns <- c("category", "subcategory", "ID", "Description", "GeneRatio", "BgRatio",
                  "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
                  "geneID", "Count")
write_result_table(kegg_df, file.path(out_root, "kegg_all.tsv"), kegg_columns)
write_result_table(kegg_sig, file.path(out_root, "kegg_significant.tsv"), kegg_columns)
kegg_status_label <- if (kegg_status == "failed") "failed" else if (nrow(kegg_df) == 0) "completed_empty" else "completed"
writeLines(c(
  paste("status", kegg_status_label),
  paste("all_terms", nrow(kegg_df)),
  paste("significant_terms", nrow(kegg_sig)),
  kegg_error
), file.path(out_root, "KEGG_STATUS.txt"))

plot_terms <- function(df, path, title, ontology_col = NULL) {
  png(path, width = 2200, height = 1600, res = 220)
  if (!nrow(df)) {
    plot.new()
    text(0.5, 0.55, "No terms passed adjusted p < 0.05", cex = 1.2)
    text(0.5, 0.45, title, cex = 1.1)
    dev.off()
    return(invisible(NULL))
  }
  df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
  if (!is.null(ontology_col) && ontology_col %in% names(df)) {
    df <- do.call(rbind, lapply(split(df, df[[ontology_col]]), head, 10))
  } else {
    df <- head(df, 15)
  }
  df$label <- factor(df$Description, levels = rev(unique(df$Description)))
  df$score <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  p <- ggplot(df, aes(x = score, y = label, size = Count, colour = if (!is.null(ontology_col) && ontology_col %in% names(df)) .data[[ontology_col]] else score)) +
    geom_point() +
    labs(x = "-log10(adjusted p-value)", y = NULL, title = title, colour = ontology_col, size = "Count") +
    theme_bw(base_size = 13)
  print(p)
  dev.off()
}

plot_terms(go_sig, file.path(out_root, "GO_recomputed.png"), "Recomputed GO enrichment", "ONTOLOGY")
plot_terms(kegg_sig, file.path(out_root, "KEGG_recomputed.png"), "Recomputed KEGG enrichment")

writeLines(c(
  "These are fresh enrichment results computed from the current experimental differential-expression tables.",
  "These are the pathway-enrichment results for the project run.",
  paste("GO terms total:", nrow(go_df)),
  paste("GO terms with adjusted p < 0.05:", nrow(go_sig)),
  paste("KEGG terms total:", nrow(kegg_df)),
  paste("KEGG terms with adjusted p < 0.05:", nrow(kegg_sig)),
  paste("KEGG status:", kegg_status_label),
  if (length(kegg_error)) paste("KEGG error:", kegg_error) else ""
), file.path(out_root, "RECOMPUTATION_NOTE.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))

cat("Recomputed enrichment written to:", out_root, "\n")
cat("GO significant terms:", nrow(go_sig), "\n")
cat("KEGG significant terms:", nrow(kegg_sig), "(status:", kegg_status_label, ")\n")
