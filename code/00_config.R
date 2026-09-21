get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (length(file_arg) == 0) {
    return(normalizePath("."))
  }
  normalizePath(sub("^--file=", "", file_arg[1]))
}

project_root <- normalizePath(file.path(dirname(get_script_path()), ".."), mustWork = TRUE)
data_root <- Sys.getenv("CAD_DATA_ROOT", unset = file.path(project_root, "data"))
data_root <- normalizePath(data_root, mustWork = TRUE)
source_root <- data_root
output_root <- Sys.getenv("CAD_OUTPUT_ROOT", unset = file.path(project_root, "generated"))
output_root <- normalizePath(output_root, mustWork = FALSE)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  gse20680_raw = file.path(data_root, "raw", "GSE20680", "geneMatrix.txt"),
  gse20680_con = file.path(data_root, "raw", "GSE20680", "s1.txt"),
  gse20680_cad = file.path(data_root, "raw", "GSE20680", "s2.txt"),
  gse20681_raw = file.path(data_root, "raw", "GSE20681", "geneMatrix.txt"),
  gse20681_con = file.path(data_root, "raw", "GSE20681", "s1.txt"),
  gse20681_cad = file.path(data_root, "raw", "GSE20681", "s2.txt"),
  gse113079_series_matrix = file.path(data_root, "raw", "GSE113079", "GSE113079_series_matrix.txt"),
  gse113079_platform = file.path(data_root, "raw", "GSE113079", "GPL20115-26806.txt"),
  batch_all = file.path(data_root, "intermediate", "batch_diff", "Batch.all.txt"),
  batch_diff = file.path(data_root, "intermediate", "batch_diff", "Batch.diff.txt"),
  diff_gene_exp = file.path(data_root, "intermediate", "batch_diff", "diffGeneExp.txt"),
  merged_before_combat = file.path(data_root, "intermediate", "normalization", "merged.before_combat.txt"),
  merged_after_combat = file.path(data_root, "intermediate", "normalization", "merged.after_combat.txt"),
  metadata = file.path(data_root, "intermediate", "normalization", "metadata.txt"),
  lasso_matrix = file.path(data_root, "intermediate", "batch_diff", "diffGeneExp.txt"),
  lasso_genes = file.path(data_root, "intermediate", "lasso", "LASSO.gene.txt"),
  rf_matrix = file.path(data_root, "intermediate", "random_forest", "diffGeneExp.txt"),
  rf_genes = file.path(data_root, "intermediate", "random_forest", "rfGenes.txt"),
  rf_importance = file.path(data_root, "intermediate", "random_forest", "geneImportance.txt"),
  svm_matrix = file.path(data_root, "intermediate", "svm", "diffGeneExp.txt"),
  svm_genes = file.path(data_root, "intermediate", "svm", "SVM-RFE.gene.txt"),
  signature_input = file.path(data_root, "intermediate", "signature", "input.txt"),
  roc_matrix = file.path(data_root, "intermediate", "roc", "diffGeneExp.txt"),
  roc_genes = file.path(data_root, "intermediate", "roc", "interGenes.txt"),
  go_bp_sig = file.path(data_root, "intermediate", "enrichment", "go.BP.sig.tsv"),
  go_cc_sig = file.path(data_root, "intermediate", "enrichment", "go.CC.sig.tsv"),
  go_mf_sig = file.path(data_root, "intermediate", "enrichment", "go.MF.sig.tsv"),
  pathway_sig = file.path(data_root, "intermediate", "enrichment", "pathway.sig.tsv")
)

required_cran_packages <- c(
  "ggplot2",
  "pheatmap",
  "glmnet",
  "randomForest",
  "e1071",
  "kernlab",
  "caret",
  "pROC",
  "VennDiagram",
  "reshape2",
  "futile.logger"
)

required_bioc_packages <- c(
  "limma",
  "sva",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

make_out_dir <- function(name) {
  path <- file.path(output_root, name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

read_expression <- function(path) {
  read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
}

read_gene_list <- function(path) {
  x <- readLines(path, warn = FALSE)
  x[nzchar(x)]
}

describe_matrix <- function(path) {
  mat <- read_expression(path)
  sprintf("%s: %d genes x %d samples", basename(path), nrow(mat), ncol(mat))
}

sample_group <- function(samples) {
  ifelse(grepl("(_con$|_Control$|Control$)", samples, ignore.case = TRUE), "Control", "CAD")
}

clean_sample_id <- function(samples) {
  sub("_(con|Control|Treat)$", "", samples, ignore.case = TRUE)
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(value) %in% c("1", "true", "yes", "y", "on")
}
