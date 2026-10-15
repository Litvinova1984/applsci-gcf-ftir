# Exploratory Chemometric Analysis of Pediatric Gingival Crevicular Fluid FTIR Measurements

**Partition Stability and Sensitivity to Sex-Associated and Clean-Strip-Aligned Variation**

Reproducibility repository for the exploratory *Applied Sciences* resubmission (15 October 2026). This project documents multivariate structure in pediatric gingival crevicular fluid (GCF) FTIR measurements, internal stability of an exploratory hierarchical partition, sex-associated spectral sensitivity, and sensitivity to an independently defined clean-strip spectral direction.

**Interpretation:** The three-group partition is a descriptive representation of the dataset. The study does **not** establish discrete biological/clinical classes or assign participant-level diagnoses based on the partition. Clean-strip-direction residualization is a sensitivity analysis, **not** a correction that recovers a purely biological spectrum.

## Analytical scope

- Active cohort: **25** pediatric participants; **12 female, 13 male**.
- Active matrix: **25 participant-level spectra × 793 variables**.
- Spectral windows: **3400–2800 and 1800–870 cm⁻¹**.
- Preprocessing: spectrum-wise standard normal variate (SNV); PCA with feature-wise mean centering and no additional feature-wise scaling.
- Baseline clustering: first **3 PCs** (**89.28%** variance), Euclidean distance, Ward.D2, descriptive **k = 3** partition (**15/4/6** participants).
- Robustness: PC-count sensitivity, leave-one-out, repeated 80% subsampling, participant bootstrap and cluster-wise Jaccard.
- Sensitivity: sex-associated variation and clean-strip-aligned direction; technical reproducibility assessed using four unused absorbent strips from one production batch.
- Random seed: **4489628**.

## Repository layout

```text
R/                         00_utils.R; 01_import_preprocess.R; 02–08 analysis scripts
analysis/                  Quarto source and rendered HTML reproducibility report
data/pediatric/             De-identified participant metadata and processed data
 data/pediatric/processed/ Frozen processed inputs and fitted analysis objects
 data/technical_controls/  Four original clean-strip DPT files
figures/                   Six main figures and four supplementary figures (PNG/PDF), plus diagnostics
supplement/                Seven supplementary CSV tables (S1–S7)
tables/                    Main tables, analysis exports, manifests, session information
tests/                     Seven R validation scripts
MANIFEST_SHA256.csv         Checksums for the archived reproducibility snapshot
README.txt                  Source snapshot note from the final package
renv.lock / renv/           R environment files retained from the previous repository version
```

## Data availability and reproducibility boundary

The repository includes processed de-identified pediatric spectra and participant-level metadata, frozen analysis objects, clean-strip technical-reference spectra, and numerical outputs. Authorization for open publication of de-identified individual-level data was confirmed by the authors. Source-linked clinical records and direct personal identifiers are not included.

**Raw OPUS.zip is not public.** Consequently, `R/01_import_preprocess.R` is retained for analytical provenance but **cannot be executed from this repository alone**. Full raw-to-processed reproduction is not claimed. The supplied frozen processed RDS/CSV inputs allow the later analysis scripts to be inspected and rerun independently of raw OPUS files, subject to the relevant R package environment.

The archive includes source inventory and import-status files with study-local, non-identifying spectral filenames. Those files are part of the documented provenance, not clinical records.

## Using the released results

The final tables, figures, and rendered report are already present. View the self-contained report at:

`analysis/revision_report.html`

Publication outputs:

- Main figures **Figure 1–Figure 6** and supplementary figures **Figure S1–Figure S4** (`figures/`, PNG/PDF).
- Three main tables (`tables/Table1_...`, `Table2_...`, `Table3_...`).
- **Seven** supplementary machine-readable tables in `supplement/`.
- Table S7: `supplement/TableS7_PCA_eigenvector_loadings_PC1_PC3.csv` (**793 rows**, five columns in the order `wavenumber, variable, PC1, PC2, PC3`).

## Re-running analyses from supplied processed inputs

Open the repository root in R/RStudio. R and the required packages are needed; Quarto must be installed separately and available on `PATH` to regenerate the HTML report. Consult the `sessionInfo` exports in `tables/` and `data/pediatric/processed/` for the software actually used in the frozen run. The existing `renv.lock` comes from the earlier repository and was retained rather than regenerated for this package; it is **not independently verified as a new lockfile for this resubmission**. Restoring or updating the environment should be done in a separate testing copy and is not required merely to view the released outputs.

For a **fresh rebuild from the provided processed inputs**, use this order (do not start at script 01 unless you possess the restricted raw OPUS archive):

```r
source("R/02_pca_clustering.R")
source("R/03_cluster_stability.R")
source("R/04_sex_sensitivity.R")
source("R/05_strip_reproducibility.R")
source("R/06_assignment_transitions.R")
source("R/07_tables_figures.R")
source("tests/validate_tables_figures_outputs.R")
source("R/08_render_reproducibility_report.R")
source("tests/validate_reproducibility_report.R")
```

This sequence **was not rerun by the final package assembler**. The author confirmed successful local execution of R07/R08 and their validators, and the delivered output files were independently checked. A complete independent R rerun from processed inputs was **not verified** as part of the package QA. Do not interpret the included frozen output files as evidence of a new end-to-end run.

## Integrity and updating

`MANIFEST_SHA256.csv` and the R07/R08 manifests document the archived reproducibility snapshot delivered for the resubmission. The archived scripts, numerical files, PNG/PDF figures, and Quarto HTML were checked against those manifests. Git repository metadata and this updated GitHub README/checklist are additional repository-level files, not entries in the archived snapshot manifest.

If any computational outputs are regenerated, rerun their validators and **update the associated manifests** before declaring a new frozen release. PDF binary checksums can change with embedded creation timestamps even when figure content is unchanged.

## Excluded analyses

The earlier public revision contained adult-reference mapping, supervised post-clustering spectral domains (D1–D5), and related outputs. They are **not part of this exploratory resubmission** and have been removed from the active release files. Historical Git commits remain available as version history; they are not part of the present frozen results.

## Citation and version

Please cite the associated *Applied Sciences* manuscript when it becomes available, and cite the **exact GitHub Release tag** used for any reproduction. The permanent release URL should be added to the manuscript's Data Availability Statement after publication of the release.
