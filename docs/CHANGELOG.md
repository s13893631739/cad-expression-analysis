# Repository preparation notes

- Copied research scripts and selected historical outputs without changing the original project.
- Changed generated output paths from `results/` to `generated/` to keep historical snapshots separate.
- Made the environment check stop with an explicit message when inputs are missing.
- Corrected transferred logistic-regression prediction to use the fitted feature scale; the corrected analysis has not been rerun.
- Excluded the old transferred-GLM performance figure and Word report to avoid mixing old conclusions with corrected code.
- Added an English README, a Chinese overview, English supporting documentation, and a base-R check of the published tables.
- Kept raw and full intermediate expression matrices outside the upload package. No publication, DOI, license, historical commit record, or full-reproducibility claim has been added.
