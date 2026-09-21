# Data availability

Source accessions: GSE20680, GSE20681, GSE113079, and platform GPL20115. These identifiers can be searched in NCBI GEO. Sample counts in the project README are project inclusion counts.

The repository contains 28 project inputs through Git LFS: 8 raw files and 20 intermediate files. `input_manifest.tsv` records their relative paths, byte sizes, and SHA-256 checksums. After cloning, run the following from the repository root:

```sh
git lfs install
git lfs pull
Rscript scripts/check_inputs.R
```

The check is read-only. It confirms that every manifest entry exists, is a downloaded file rather than a Git LFS pointer, and matches its recorded size and checksum. `git lfs ls-files` should list the same 28 paths.

## Local layout

The expected `data/raw/` and `data/intermediate/` directories are included in the repository through Git LFS. Preserve file names and subfolders when replacing an input. Do not commit alternate inputs without updating the manifest and documenting the change.

The raw directory contains the discovery gene matrices, group lists, external series matrix, and platform annotation. The intermediate directory contains normalization matrices, sample metadata, differential-expression tables, model input matrices, selected gene lists, and enrichment tables. Exact expected paths are listed in [the manifest](input_manifest.tsv).

## Reproduction from public records

The repository includes the restored project-specific input files, but it does not include every historical step used to obtain `geneMatrix.txt` and the intermediate inputs from GEO. The full experimental script can operate on the included inputs; a clean GEO-download-to-input reconstruction would still require the original acquisition, annotation, grouping, and preprocessing steps.

The small published tables in `results/` are sufficient for `scripts/check_results.R`, not for retraining the models.
