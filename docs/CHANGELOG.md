# Repository preparation notes

- Copied research scripts and selected historical outputs without changing the original project.
- Changed generated output paths from `results/` to `generated/` to keep historical snapshots separate.
- Made the environment check stop with an explicit message when inputs are missing.
- Corrected transferred logistic-regression prediction and reran the saved-input evaluation with training-cohort direction and scaling.
- Excluded the old transferred-GLM performance figure and Word report to avoid mixing old conclusions with corrected code.
- Added an English README, a Chinese overview, English supporting documentation, and a base-R check of the published tables.
- The initial repository snapshot kept raw and full intermediate expression matrices outside the upload package. No publication, DOI, license, or historical commit record has been added.
- The restored-input experimental run is recorded separately because its current-environment counts differ from the historical snapshot; historical `results/` files remain unchanged.
- Restored the 28 manifest inputs with Git LFS and added a separate `run_full_experiment.R` entry point.
- Made the saved-input runner record `sessionInfo()` and run metadata, and made external ROC direction/scaling use training-cohort values.
