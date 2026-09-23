source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1]))), "00_config.R"))

cat("R version:", R.version.string, "\n")
cat("Data:", data_root, "\n")
cat("Output:", output_root, "\n\n")

required_inputs <- c(
  "gse20680_raw", "gse20680_con", "gse20680_cad",
  "gse20681_raw", "gse20681_con", "gse20681_cad",
  "gse113079_series_matrix", "gse113079_platform"
)
missing_files <- required_inputs[!file.exists(unlist(paths[required_inputs]))]
if (length(missing_files) == 0) {
  cat("All known input files exist.\n\n")
} else {
  cat("Missing input files:\n")
  for (name in missing_files) {
    cat(" -", name, ":", paths[[name]], "\n")
  }
  stop("Required local inputs are missing. See data/README.md and data/input_manifest.tsv.")
}

is_lfs_pointer <- function(path) {
  if (!file.exists(path)) return(FALSE)
  first_line <- readLines(path, n = 1, warn = FALSE)
  length(first_line) == 1 && identical(first_line, "version https://git-lfs.github.com/spec/v1")
}

lfs_pointers <- required_inputs[vapply(paths[required_inputs], is_lfs_pointer, logical(1))]
if (length(lfs_pointers) > 0) {
  stop(
    "Git LFS objects are not downloaded for: ",
    paste(lfs_pointers, collapse = ", "),
    ". Run `git lfs pull` and retry."
  )
}

for (name in c("gse20680_raw", "gse20681_raw")) {
  cat(describe_matrix(paths[[name]]), "\n")
}
cat("\n")

all_packages <- c(required_cran_packages, required_bioc_packages)
installed <- rownames(installed.packages())
cat("R package status:\n")
for (pkg in all_packages) {
  status <- if (pkg %in% installed) "OK" else "MISSING"
  cat(sprintf("%-7s %s\n", status, pkg))
}
