# Running the CAD expression-analysis project

Run all commands from the repository root. The canonical workflow is `code/run_project.R`; it is the single end-to-end project entry point.

## 1. Hydrate and verify inputs

```sh
git lfs install
git lfs pull
Rscript scripts/check_inputs.R
```

The input checker verifies the files listed in `data/input_manifest.tsv`, including byte sizes and SHA-256 checksums. A missing file or Git LFS pointer stops the run.

## 2. Run the complete project

```sh
Rscript code/run_project.R
```

The runner performs, in order:

1. environment and input checks;
2. discovery normalization, ComBat batch correction, QC and limma differential expression;
3. descriptive LASSO, random-forest and SVM-RFE candidate screens;
4. fresh GO and KEGG enrichment from the current differential-expression tables;
5. nested Elastic Net model assessment and final model fitting;
6. one frozen-model exploratory evaluation in GSE113079;
7. a unified run record and project summary.

The main output tree is `generated/project/`. The repository-level display files are `results/project_summary.tsv`, `results/project_summary.md`, and `results/predictive_signature.tsv`.

## 3. Main outputs

- `generated/project/20_full_experiment/`: normalization, QC, differential expression and descriptive candidate screens.
- `generated/project/12_enrichment_recomputed/`: current GO/KEGG results and `KEGG_STATUS.txt`.
- `generated/project/50_nested_cv_signature/`: outer predictions, metrics, gene frequencies, coefficients and preprocessing parameters.
- `generated/project/51_nested_signature_external/`: frozen-model exploratory predictions, probe mapping and external metrics.
- `generated/project/RUN_METADATA.txt`: canonical run description.
- `generated/project/sessionInfo.txt`: R and package session information.

The external evaluation records any final model genes that lack GPL20115 probes in `unmapped_frozen_genes.txt` and omits them from that exploratory score. The discovery model itself remains fully recorded in `predictive_signature.tsv`.

## 4. Reproducibility notes

The predictive stage uses the fresh normalization artifacts created in the same run through `CAD_DISCOVERY_ROOT`. Feature screening and training-fold scaling are fitted inside the nested model-development splits. The supplied normalized discovery matrix is retained as the common project input because ComBat does not provide a general frozen prediction transform for an unseen cohort.

KEGG depends on access to the KEGG REST service. If the service is unavailable, the run records `status failed` and preserves the GO results.
