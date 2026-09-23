# Apply the frozen nested-CV Elastic Net signature to GSE113079.
#
# This script is exploratory because GSE113079 is an external cohort. It never
# uses its labels to choose genes,
# scaling, direction, alpha, or lambda. The discovery model and preprocessing
# parameters are read from generated/50_nested_cv_signature.

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- command_args[grepl("^--file=", command_args)]
if (length(file_arg) != 1) stop("Run this script with Rscript code/51_apply_nested_signature_external.R")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
source(file.path(script_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(pROC)
})

model_root <- file.path(output_root, "50_nested_cv_signature")
if (!file.exists(file.path(model_root, "final_coefficients.tsv"))) {
  stop("Run code/50_nested_cv_signature.R first.")
}
out_root <- file.path(output_root, "51_nested_signature_external")
if (dir.exists(out_root)) unlink(list.files(out_root, full.names = TRUE, all.files = TRUE, no.. = TRUE), recursive = TRUE, force = TRUE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

coefficients <- read.table(file.path(model_root, "final_coefficients.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
coefficients <- coefficients[coefficients$Gene != "(Intercept)" & coefficients$coefficient != 0, , drop = FALSE]
intercept_path <- file.path(model_root, "final_coefficients.tsv")
# The intercept is retained in the TSV only when nonzero; read it from the
# same frozen fit artifact and fail explicitly if it was omitted.
all_coefficients <- read.table(intercept_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
intercept <- all_coefficients$coefficient[match("(Intercept)", all_coefficients$Gene)]
if (!length(intercept) || !is.finite(intercept)) stop("Frozen coefficient file does not contain an intercept.")
scaler <- read.table(file.path(model_root, "final_preprocessing_params.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
scaler <- scaler[scaler$Gene %in% coefficients$Gene, , drop = FALSE]
if (!nrow(scaler)) stop("No frozen nonzero genes have external probe mappings.")

split_gene_symbols <- function(x) {
  if (is.na(x) || !nzchar(x) || x %in% c("NA", "None")) return(character(0))
  trimws(unlist(strsplit(x, "[,;|/]+")))
}
external <- readLines(paths$gse113079_series_matrix, warn = FALSE)
begin <- grep("^!series_matrix_table_begin", external)
end <- grep("^!series_matrix_table_end", external)
expr_df <- read.table(paths$gse113079_series_matrix, skip = begin, nrows = end - begin - 2, header = TRUE, sep = "\t", quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_df) <- expr_df[[1]]
expr <- data.matrix(expr_df[, -1, drop = FALSE])
title_idx <- grep("^!Sample_title", external)
sample_titles <- gsub("^\"|\"$", "", strsplit(external[title_idx], "\t", fixed = TRUE)[[1]][-1])
group <- ifelse(grepl("healthy control", sample_titles, ignore.case = TRUE), "Control", "CAD")
if (any(!grepl("healthy control|CAD", sample_titles, ignore.case = TRUE))) stop("Unrecognized external sample title.")
group <- factor(group, levels = c("Control", "CAD"))
platform <- read.delim(paths$gse113079_platform, skip = 41, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

mapped <- lapply(scaler$Gene, function(gene) {
  hits <- vapply(platform$GeneSymbol, function(x) gene %in% split_gene_symbols(x), logical(1))
  probes <- intersect(unique(na.omit(platform$ID[hits])), rownames(expr))
  if (!length(probes)) return(list(Gene = gene, probes = character(), values = rep(NA_real_, ncol(expr))))
  values <- if (length(probes) == 1) as.numeric(expr[probes, ]) else colMeans(expr[probes, , drop = FALSE], na.rm = TRUE)
  list(Gene = gene, probes = probes, values = values)
})
missing_genes <- vapply(mapped, function(x) length(x$probes) == 0, logical(1))
missing_names <- vapply(mapped[missing_genes], `[[`, character(1), "Gene")
if (any(missing_genes)) {
  warning("Frozen model genes without GPL20115 probes were omitted from the exploratory external score: ", paste(missing_names, collapse = ";"))
  writeLines(missing_names, file.path(out_root, "unmapped_frozen_genes.txt"))
}
mapped <- mapped[!missing_genes]
if (!length(mapped)) stop("None of the frozen nonzero genes map to GSE113079 probes.")
external_expr <- do.call(cbind, lapply(mapped, `[[`, "values"))
colnames(external_expr) <- vapply(mapped, `[[`, character(1), "Gene")
scaler <- scaler[match(colnames(external_expr), scaler$Gene), , drop = FALSE]
coefficients <- coefficients[match(colnames(external_expr), coefficients$Gene), , drop = FALSE]
probe_map <- data.frame(Gene = colnames(external_expr), Probes = vapply(mapped, function(x) paste(x$probes, collapse = ";"), character(1)), nProbes = vapply(mapped, function(x) length(x$probes), integer(1)), stringsAsFactors = FALSE)
external_scaled <- sweep(sweep(external_expr, 2, scaler$center, "-"), 2, scaler$scale, "/")
coef_vector <- setNames(coefficients$coefficient, coefficients$Gene)
score <- as.numeric(intercept + external_scaled[, names(coef_vector), drop = FALSE] %*% coef_vector)
probability <- plogis(score)
roc_obj <- roc(group, probability, levels = c("Control", "CAD"), direction = "<", quiet = TRUE)
auc_value <- as.numeric(auc(roc_obj))

precision_recall <- function(response, score) {
  positive <- response == "CAD"
  order_index <- order(score, decreasing = TRUE)
  truth <- positive[order_index]
  tp <- cumsum(truth); fp <- cumsum(!truth)
  recall <- tp / sum(truth); precision <- tp / (tp + fp)
  sum(diff(c(0, recall)) * precision)
}
summary <- data.frame(item = c("dataset", "samples_total", "samples_cad", "samples_control", "frozen_gene_count", "mapped_frozen_gene_count", "unmapped_frozen_genes", "frozen_alpha", "external_roc_auc", "external_pr_auc", "model_source", "label_usage"), value = c("GSE113079", length(group), sum(group == "CAD"), sum(group == "Control"), nrow(all_coefficients[all_coefficients$Gene != "(Intercept)" & all_coefficients$coefficient != 0, , drop = FALSE]), nrow(coefficients), paste(setdiff(all_coefficients$Gene[all_coefficients$Gene != "(Intercept)" & all_coefficients$coefficient != 0], coefficients$Gene), collapse = ";"), read.table(file.path(model_root, "summary.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)$value[3], sprintf("%.4f", auc_value), sprintf("%.4f", precision_recall(group, probability)), "generated/50_nested_cv_signature", "labels used only for final exploratory metric"), stringsAsFactors = FALSE)
write.table(data.frame(sample = colnames(expr), group = group, score = score, probability = probability), file.path(out_root, "external_predictions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(probe_map, file.path(out_root, "probe_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(out_root, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
writeLines(c("Frozen model from nested discovery-cohort CV.", "No GSE113079 labels were used for model selection or preprocessing.", "This is an exploratory external evaluation, not a blinded confirmatory test.", paste("Mapped frozen genes used:", nrow(coefficients)), paste("Unmapped frozen genes omitted:", length(missing_names)), paste("ROC-AUC:", sprintf("%.4f", auc_value)), paste("PR-AUC:", sprintf("%.4f", precision_recall(group, probability)))), file.path(out_root, "RUN_METADATA.txt"))
cat(sprintf("Frozen nested signature external evaluation written to %s\nROC-AUC: %.4f; PR-AUC: %.4f\n", out_root, auc_value, precision_recall(group, probability)))
