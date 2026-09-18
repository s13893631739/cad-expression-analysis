source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1]))), "00_config.R"))

cat("R version:", R.version.string, "\n")
cat("Data:", data_root, "\n")
cat("Output:", output_root, "\n\n")

missing_files <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_files) == 0) {
  cat("All known input files exist.\n\n")
} else {
  cat("Missing input files:\n")
  for (name in missing_files) {
    cat(" -", name, ":", paths[[name]], "\n")
  }
  stop("Required local inputs are missing. See data/README.md and data/input_manifest.tsv.")
}

for (name in c("gse20680_raw", "gse20681_raw", "diff_gene_exp", "rf_matrix", "roc_matrix")) {
  cat(describe_matrix(paths[[name]]), "\n")
}
cat("\n")

lasso_genes <- read_gene_list(paths$lasso_genes)
rf_genes <- read_gene_list(paths$rf_genes)
svm_genes <- read_gene_list(paths$svm_genes)
intersect_genes <- Reduce(intersect, list(lasso_genes, rf_genes, svm_genes))

for (name in c("lasso_genes", "rf_genes", "svm_genes", "roc_genes")) {
  cat(name, ":", length(read_gene_list(paths[[name]])), "genes\n")
}
cat("intersect_genes :", length(intersect_genes), "genes\n")
cat("\n")

all_packages <- c(required_cran_packages, required_bioc_packages)
installed <- rownames(installed.packages())
cat("R package status:\n")
for (pkg in all_packages) {
  status <- if (pkg %in% installed) "OK" else "MISSING"
  cat(sprintf("%-7s %s\n", status, pkg))
}
