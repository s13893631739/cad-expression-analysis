source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1]))), "00_config.R"))

repos <- "https://cloud.r-project.org"
installed <- rownames(installed.packages())

missing_cran <- setdiff(required_cran_packages, installed)
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = repos)
}

if (!"BiocManager" %in% rownames(installed.packages())) {
  install.packages("BiocManager", repos = repos)
}

installed <- rownames(installed.packages())
missing_bioc <- setdiff(required_bioc_packages, installed)
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

cat("Package installation check complete.\n")
