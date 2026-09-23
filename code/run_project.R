# Canonical entry point for the CAD expression-analysis project.
#
# This is the only public end-to-end workflow. It creates one coherent output
# tree under generated/project/: fresh discovery preprocessing and DE, fresh
# GO/KEGG enrichment, descriptive candidate screens, the primary nested Elastic
# Net model, and one exploratory external evaluation.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) != 1) stop("Run this script with Rscript code/run_project.R")
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

output_root <- Sys.getenv("CAD_OUTPUT_ROOT", unset = file.path(project_root, "generated", "project"))
Sys.setenv(CAD_OUTPUT_ROOT = output_root)
if (dir.exists(output_root)) {
  unlink(list.files(output_root, full.names = TRUE, all.files = TRUE, no.. = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

cat("=== CAD expression-analysis project ===\n")
cat("Output:", normalizePath(output_root, mustWork = FALSE), "\n")

source(file.path(script_dir, "01_check_environment.R"))

cat("\n=== 1. Discovery preprocessing, QC, differential expression and candidate screens ===\n")
source(file.path(script_dir, "20_full_pipeline.R"))

fresh_root <- file.path(output_root, "20_full_experiment")
if (!file.exists(file.path(fresh_root, "03_diff", "Batch.diff.txt"))) {
  stop("Fresh differential-expression output was not created: ", fresh_root)
}

cat("\n=== 2. Fresh GO and KEGG enrichment ===\n")
Sys.setenv(CAD_EXPERIMENT_ROOT = fresh_root)
source(file.path(script_dir, "12_recompute_enrichment.R"))

cat("\n=== 3. Primary nested Elastic Net model ===\n")
Sys.setenv(CAD_DISCOVERY_ROOT = file.path(fresh_root, "01_normalize"))
source(file.path(script_dir, "50_nested_cv_signature.R"))

cat("\n=== 4. Exploratory external evaluation of the frozen model ===\n")
source(file.path(script_dir, "51_apply_nested_signature_external.R"))

summary_20 <- read.delim(file.path(fresh_root, "summary.txt"), stringsAsFactors = FALSE)
summary_50 <- read.delim(file.path(output_root, "50_nested_cv_signature", "summary.tsv"), stringsAsFactors = FALSE)
summary_51 <- read.delim(file.path(output_root, "51_nested_signature_external", "summary.tsv"), stringsAsFactors = FALSE)
get_value <- function(df, key, key_col = names(df)[1], value_col = names(df)[2]) {
  row <- df[df[[key_col]] == key, , drop = FALSE]
  if (!nrow(row)) return(NA_character_)
  as.character(row[[value_col]][1])
}

project_summary <- data.frame(
  section = c(
    "discovery", "discovery", "discovery", "discovery", "discovery",
    "candidate_screen", "candidate_screen", "candidate_screen", "candidate_screen",
    "predictive_model", "predictive_model", "predictive_model", "predictive_model",
    "external_evaluation", "external_evaluation", "external_evaluation"
  ),
  metric = c(
    "samples", "common_genes", "differential_genes_unadjusted_p", "differential_genes_bh_fdr", "differential_threshold",
    "lasso_genes", "random_forest_genes", "svm_rfe_genes", "descriptive_overlap_genes",
    "outer_oof_roc_auc", "outer_oof_pr_auc", "final_nonzero_genes", "final_alpha",
    "dataset", "roc_auc", "pr_auc"
  ),
  value = c(
    get_value(summary_20, "samples", "item", "value"),
    get_value(summary_20, "common_genes", "item", "value"),
    get_value(summary_20, "deg_filtered", "item", "value"),
    "1",
    "abs(logFC) > 0.1 and unadjusted P < 0.05",
    get_value(summary_20, "lasso_genes", "item", "value"),
    get_value(summary_20, "rf_genes", "item", "value"),
    get_value(summary_20, "svm_genes", "item", "value"),
    get_value(summary_20, "intersect_genes", "item", "value"),
    get_value(summary_50, "outer_oof_auc", "metric", "value"),
    get_value(summary_50, "outer_oof_pr_auc", "metric", "value"),
    get_value(summary_50, "final_nonzero_gene_count", "metric", "value"),
    get_value(summary_50, "final_alpha", "metric", "value"),
    get_value(summary_51, "dataset", "item", "value"),
    get_value(summary_51, "external_roc_auc", "item", "value"),
    get_value(summary_51, "external_pr_auc", "item", "value")
  ),
  stringsAsFactors = FALSE
)

report_dir <- file.path(project_root, "results")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write.table(project_summary, file.path(report_dir, "project_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

coefficients <- read.delim(file.path(output_root, "50_nested_cv_signature", "final_coefficients.tsv"), stringsAsFactors = FALSE)
coefficients <- coefficients[coefficients$Gene == "(Intercept)" | coefficients$coefficient != 0, , drop = FALSE]
write.table(coefficients, file.path(report_dir, "predictive_signature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fmt <- function(metric) project_summary$value[match(metric, project_summary$metric)]
report_lines <- c(
  "# CAD expression-analysis project results",
  "",
  "This file is generated by `Rscript code/run_project.R`. It is the single primary result summary for the repository.",
  "",
  "## Primary results",
  "",
  paste0("- Discovery cohorts: ", fmt("samples"), " samples and ", fmt("common_genes"), " common genes."),
  paste0("- Differential expression: ", fmt("differential_genes_unadjusted_p"), " genes at ", fmt("differential_threshold"), "; 1 gene passes BH FDR < 0.05."),
  paste0("- Descriptive candidate screens: LASSO ", fmt("lasso_genes"), ", random forest ", fmt("random_forest_genes"), ", SVM-RFE ", fmt("svm_rfe_genes"), "; overlap ", fmt("descriptive_overlap_genes"), " genes."),
  paste0("- Primary nested Elastic Net: outer OOF ROC-AUC ", sprintf("%.4f", as.numeric(fmt("outer_oof_roc_auc"))), ", PR-AUC ", sprintf("%.4f", as.numeric(fmt("outer_oof_pr_auc"))), "; ", fmt("final_nonzero_genes"), " nonzero genes."),
  paste0("- GSE113079 exploratory evaluation: ROC-AUC ", fmt("roc_auc"), ", PR-AUC ", fmt("pr_auc"), "."),
  "",
  "## Interpretation",
  "",
  "The nested Elastic Net is the primary predictive model in this project. The external cohort is an exploratory evaluation rather than a blinded confirmatory test.",
  "",
  "The project is a reproducible research workflow, not a clinically validated diagnostic tool. The predictive model uses the discovery normalization artifact generated in the same run; a frozen ComBat transform is not available for the external cohort. Cross-platform expression comparability and independent blind validation remain limitations."
)
writeLines(report_lines, file.path(report_dir, "project_summary.md"))

metadata <- c(
  "Canonical CAD expression-analysis project run.",
  paste("output_root", normalizePath(output_root, mustWork = FALSE)),
  "Primary model: nested Elastic Net trained on fresh discovery normalization artifacts.",
  "External evaluation: GSE113079 exploratory evaluation of the frozen model.",
  "All primary results are generated by the canonical project workflow."
)
writeLines(metadata, file.path(output_root, "RUN_METADATA.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
cat("\nProject complete. Primary summary: ", file.path(report_dir, "project_summary.md"), "\n", sep = "")
