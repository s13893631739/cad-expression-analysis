# CAD Transcriptomics and Predictive Signature

An end-to-end, reproducible R project for coronary artery disease gene-expression analysis, pathway enrichment and predictive signature modeling.

## Project overview

The workflow integrates two discovery cohorts and one external cohort to characterize CAD-associated transcriptional changes and build an exploratory gene-expression signature. It runs from input validation through quality control, normalization, batch correction, differential expression, pathway enrichment, feature screening, nested Elastic Net modeling and frozen external evaluation.

## Key results

- **Discovery data:** 337 samples and 19,261 shared genes from GSE20680 and GSE20681.
- **Differential expression:** 286 genes at `abs(logFC) > 0.1` and nominal `P < 0.05`; one gene passes BH FDR < 0.05.
- **Candidate screening:** LASSO 32 genes, random forest 39 genes, SVM-RFE 31 genes, with a 15-gene descriptive overlap.
- **Predictive model:** nested Elastic Net outer OOF ROC-AUC **0.5505**, PR-AUC **0.5969**, with 6 nonzero genes in the final model.
- **External evaluation:** frozen-model exploratory evaluation on GSE113079, ROC-AUC **0.4415** and PR-AUC **0.6481**.

Main outputs:

- [Project summary](results/project_summary.md)
- [Metrics table](results/project_summary.tsv)
- [Predictive signature coefficients](results/predictive_signature.tsv)
- [Methods](docs/METHODS.md)
- [Running instructions](docs/RUNNING.md)

Figures:

- [Batch-corrected PCA](figures/pca_after_batch.png)
- [Differential-expression volcano plot](figures/volcano.png)
- [GO enrichment](figures/go_enrichment.png)
- [Nested-model ROC](figures/nested_model_roc.png)
- [External-evaluation ROC](figures/external_model_roc.png)

## Data

| Dataset | Role | CAD | Control |
|---|---|---:|---:|
| GSE20680 | Discovery | 87 | 52 |
| GSE20681 | Discovery | 99 | 99 |
| GSE113079 | External evaluation | 93 | 48 |

The project uses raw GEO matrices and sample-group files stored with Git LFS. GPL20115 provides the external probe annotation.

## Analysis workflow

1. Verify raw inputs and SHA-256 checksums.
2. Normalize the discovery matrices, merge shared genes and correct dataset batch effects with ComBat.
3. Generate QC plots and run limma CAD-versus-control differential analysis.
4. Recompute GO and KEGG enrichment with an explicit gene universe and multiple-testing correction.
5. Generate descriptive LASSO, random-forest and SVM-RFE candidate screens.
6. Fit nested Elastic Net models; feature filtering, centering/scaling and alpha selection are performed inside training folds.
7. Freeze the final discovery model and evaluate it on GSE113079 with fixed coefficients and probe mapping.

## Technical stack

R, limma, sva/ComBat, clusterProfiler, org.Hs.eg.db, glmnet, randomForest, caret, e1071, ggplot2, Git LFS and GitHub Actions.

## Reproduce the project

```sh
git lfs install
git lfs pull
Rscript scripts/check_inputs.R
Rscript code/run_project.R
```

The single public entry point writes the full run to `generated/project/` and records `sessionInfo()`. Compact results are committed under `results/` so the project can be reviewed directly on GitHub.

## Repository structure

```text
code/run_project.R       end-to-end project entry point
code/                    analysis stages
scripts/check_inputs.R   input verification
data/raw/                GEO matrices, sample groups and platform annotation
results/                 project metrics and predictive signature
figures/                 project figure notes
docs/                    methods and reproducibility notes
generated/project/       generated stage outputs
```

## Limitations

This is an exploratory research workflow, not a clinically validated diagnostic tool. The external cohort is an exploratory evaluation rather than a blinded confirmatory cohort. Cross-platform expression comparability, calibrated probabilities, confidence intervals and a completely untouched confirmation cohort remain open validation items.

## Resume-ready description

Built a reproducible R pipeline for CAD transcriptomics using 337 samples from two GEO cohorts; implemented normalization, ComBat batch correction, limma differential expression, GO/KEGG enrichment and multi-method candidate screening; developed a training-fold-controlled nested Elastic Net signature with OOF ROC-AUC 0.5505 and PR-AUC 0.5969; packaged the workflow with checksum-verified inputs, reports and one-command execution.

[中文说明](README.zh-CN.md)
