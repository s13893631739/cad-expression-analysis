# Methods and interpretation

## Study design and preprocessing

The project combines GSE20680 and GSE20681 as discovery cohorts (337 samples, 19,261 shared genes). Each discovery matrix is normalized with `normalizeBetweenArrays`, merged on shared gene identifiers and adjusted with ComBat using dataset as the batch variable and disease group in the model matrix. PCA and boxplots are generated before and after batch correction.

The predictive stage uses the normalized discovery artifact created by the same project run. The raw GEO inputs, package versions and session information are recorded so the analysis can be rerun from the repository.

## Differential expression and enrichment

A limma contrast compares CAD/Treat with Control. The descriptive differential-expression screen uses `abs(logFC) > 0.1` and unadjusted `P < 0.05`. The current run contains 286 genes under this threshold; one passes BH-adjusted FDR < 0.05. The threshold is reported explicitly and should not be described as an FDR-significant gene set.

GO and KEGG are computed from the project differential-expression table using `org.Hs.eg.db` and clusterProfiler, a SYMBOL-to-ENTREZID mapping, gene-set sizes 10–500 and BH-adjusted `p < 0.05` for significance. KEGG service failures are recorded in `KEGG_STATUS.txt`.

## Candidate gene screens

LASSO, random forest and SVM-RFE are descriptive candidate screens. The current run produces 32, 39 and 31 genes, respectively, with a 15-gene overlap. These lists describe the exploratory feature-selection stage; their overlap is not the final predictive signature.

## Primary predictive model

The primary predictive model is an Elastic Net logistic model assessed with nested cross-validation in the discovery cohorts. Candidate genes are restricted to genes with a fixed mapping in the GPL20115 annotation. Within each training split, genes are screened, centering and scaling parameters are fitted, and the Elastic Net alpha is selected by inner validation. The outer predictions provide the primary internal estimate. The final model is refit after model assessment and its coefficients and preprocessing parameters are saved.

The current outer out-of-fold results are ROC-AUC 0.5505 and PR-AUC 0.5969, with six nonzero final model genes. These are model-development estimates within the discovery cohorts, not proof of clinical performance.

## External evaluation

The final discovery model is applied once to GSE113079 using the frozen coefficients and discovery preprocessing parameters. Probe-to-gene mapping uses GPL20115. Genes without a matching expression probe are recorded and omitted from the exploratory external score. The current exploratory result is ROC-AUC 0.4415 and PR-AUC 0.6481.

This evaluation is exploratory rather than a blinded confirmatory test. Cross-platform expression comparability remains unresolved, and the external result must not be used to tune the model.

## Interpretation and scope

This is a reproducible exploratory research workflow, not a clinically validated diagnostic tool. A completely untouched confirmation cohort, calibrated probability estimates, uncertainty intervals and cross-study leave-one-study-out validation remain future work.
