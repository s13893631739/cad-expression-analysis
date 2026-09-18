# Methods and interpretation

## Preprocessing and differential expression

The discovery datasets are normalized separately with `normalizeBetweenArrays`, combined by shared gene identifiers, and adjusted with ComBat using dataset as the batch variable and group in the model matrix. Boxplots and PCA are used for inspection.

The limma contrast compares CAD/Treat with Control. The saved candidate set contains 284 genes (221 up, 63 down), selected with absolute logFC > 0.1 and unadjusted P < 0.05. Only one gene in this set meets BH-adjusted P < 0.05. The candidate set should not be described as 284 FDR-significant genes.

GO and KEGG figures use existing result tables. Complete upstream enrichment calculations, including the background and all parameters, are not provided.

## Feature selection

The saved lists contain 32 LASSO, 39 random forest, and 31 SVM-RFE genes, with an intersection of 10. The code includes binomial LASSO with 10-fold cross-validation, a 500-tree random forest, and radial-kernel SVM-RFE with 10-fold cross-validation.

Selection is not fully independent across methods. The experimental path feeds LASSO candidates into SVM-RFE and retains fixed target counts for random forest and SVM. The figure-generation path consumes historical feature lists. Consequently, the published overlap describes saved selections rather than an independently validated optimal signature.

The saved SVM table reports Accuracy 0.691021 at 31 variables. Screening before cross-validation can bias this estimate; it is not independent test accuracy.

## External evaluation

GSE113079 includes 93 CAD and 48 control samples in this analysis. Probes are mapped using GPL20115 GeneSymbol annotations, with mean aggregation for genes having multiple probes. MAPK8IP1 has a saved AUC of 0.765681 and CLEC4D 0.683020.

ROC direction is chosen automatically on the evaluated data, so AUC alone does not establish consistent biological direction across cohorts. The historical weighted 10-gene score has AUC 0.539 and uses cohort-level external standardization. It is not a frozen single-patient deployment procedure.

The original transferred logistic-regression score combined standardized external features with coefficients fitted on unstandardized discovery features. Repository preparation changed prediction to use the fitted model's feature scale consistently. This has not been rerun, and cross-platform comparability remains unresolved. The original transfer AUC and its figure are excluded from the published snapshots.

## Scope

This is an exploratory analysis, not a clinically validated diagnostic tool. Full reproducibility requires the project-specific raw and intermediate inputs; a complete GEO-download-to-results workflow and historical dependency versions are not available.
