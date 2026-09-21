# CAD expression analysis

[中文说明](README.zh-CN.md)

R scripts and results from an independent study of gene expression in coronary artery disease (CAD). The analysis combines two discovery datasets, screens candidate genes, and examines their performance in a separate dataset.

## Data and methods

| Dataset | Role | CAD | Control |
|---|---|---:|---:|
| GSE20680 | Discovery | 87 | 52 |
| GSE20681 | Discovery | 99 | 99 |
| GSE113079 | External evaluation | 93 | 48 |

These are the samples included in this project. The discovery matrices contain 19,261 genes across 337 samples; external probe annotation uses GPL20115.

The analysis uses expression normalization, ComBat batch correction, PCA, limma differential expression analysis, and LASSO/random forest/SVM-RFE candidate screening. Existing GO and KEGG result tables are used for enrichment interpretation and plotting. External evaluation includes probe mapping and ROC analysis. The restored raw and intermediate inputs are versioned with Git LFS.

## Selected results

The saved screening lists share 10 candidate genes. In GSE113079, MAPK8IP1 has a single-gene AUC of 0.766; the historical 10-gene weighted score has an AUC of 0.539. These are exploratory results, with limited discrimination for the combined score.

![External single-gene ROC curves in GSE113079](figures/external_single_gene_roc.png)

[Result tables and figure guide](results/README.md) · [Methods and limitations](docs/METHODS.md)

## Getting started

To inspect the included gene lists and AUC table, run this from the repository root with R installed:

```sh
Rscript scripts/check_results.R
```

This uses base R and the included files only. It checks the saved gene intersection and prints saved AUC values; it does not retrain models or recompute ROC curves.

To rerun the saved-input evaluation, install Git LFS and the listed R packages, then run `Rscript code/run_final.R`. The full experimental path is available as `Rscript code/run_full_experiment.R`; its outputs remain separate from the saved-input evaluation. See [data availability](data/README.md) and [running the analysis](docs/RUNNING.md).

## Files

```text
code/       Analysis scripts
scripts/    Check the included result tables
results/    Saved numerical results and gene lists
figures/    Selected figures from the original analysis
data/       Input requirements and file checksums
docs/       Methods and running instructions
```

## Author

Jialiang Sun (孙嘉良). Data preparation, analysis, interpretation, and report writing were carried out as an independent project.

For questions about the analysis, use this repository's Issues page when available. No software license has been selected. Third-party data and dependencies retain their respective terms.
