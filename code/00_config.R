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
  gse113079_platform = file.path(data_root, "raw", "GSE113079", "GPL20115-26806.txt")
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
