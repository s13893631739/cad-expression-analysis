# Project data

The project uses GEO accessions GSE20680, GSE20681 and GSE113079, with GPL20115 for external probe annotation. Raw matrices and sample-group files are stored with Git LFS.

```sh
git lfs install
git lfs pull
Rscript scripts/check_inputs.R
```

The checksum manifest lists the raw files required by the canonical workflow. The analysis starts from these inputs and writes all derived data under `generated/project/`.
