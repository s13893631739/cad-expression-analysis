# Data availability

Source accessions: GSE20680, GSE20681, GSE113079, and platform GPL20115. These identifiers can be searched in NCBI GEO. Sample counts in the project README are project inclusion counts.

The repository does not contain the raw expression matrices or complete intermediate data. `input_manifest.tsv` records the original local data paths, byte sizes, and SHA-256 checksums, not their contents.

## Local layout

Restore the original project's `data/raw/` and `data/intermediate/` directories here, preserving file names and subfolders. These directories are excluded by `.gitignore`. Browser uploads do not apply that exclusion automatically.

The raw directory contains the discovery gene matrices, group lists, external series matrix and platform annotation. The intermediate directory contains normalization matrices, sample metadata, differential-expression tables, model input matrices, selected gene lists and enrichment tables. Exact expected paths are listed in [the manifest](input_manifest.tsv).

## Reproduction from public records

The original package does not include all steps used to obtain `geneMatrix.txt` and the historical intermediate inputs from GEO. Downloading the source records alone is therefore insufficient to reproduce the saved results. A future reproducible release would need those acquisition, annotation, grouping and preprocessing steps, with their parameters and dependency versions.

The small published tables in `results/` are sufficient for `scripts/check_results.R`, not for retraining the models.
