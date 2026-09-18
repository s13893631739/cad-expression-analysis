# Inspect the included result tables using base R only.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script), ".."))
read_genes <- function(path) {
  x <- readLines(file.path(root, path), warn = FALSE)
  unique(x[nzchar(x)])
}
sets <- list(
  LASSO = read_genes("results/05_lasso/LASSO.gene.txt"),
  RandomForest = read_genes("results/06_random_forest/rfGenes.txt"),
  SVM_RFE = read_genes("results/07_svm/SVM-RFE.gene.txt")
)
shared <- Reduce(intersect, sets)
saved <- read_genes("results/08_venn/intersectGenes.txt")
if (!setequal(shared, saved)) stop("Saved intersection differs from the source lists.")
auc <- read.delim(file.path(root, "results/11_external_validation/roc_auc_table.txt"))
if (anyDuplicated(auc$Gene) || !setequal(auc$Gene, shared)) {
  stop("External AUC table does not match the candidate gene set.")
}
if (any(!is.finite(auc$AUC) | auc$AUC < 0 | auc$AUC > 1)) stop("Invalid AUC values.")
cat("Saved gene lists:\n")
print(vapply(sets, length, integer(1)))
cat("Shared candidates:", length(shared), "\n")
cat("External single-gene AUC (saved values, not recomputed):\n")
print(auc[order(-auc$AUC), c("Gene", "AUC")], row.names = FALSE)
cat("Result-table consistency checks passed. No model was trained.\n")
