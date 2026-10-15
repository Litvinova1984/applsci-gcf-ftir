Applied Sciences reproducibility snapshot — 15 October 2026
Status: READY FOR GITHUB RELEASE (technical package assembly)

Manuscript title:
Exploratory Chemometric Analysis of Pediatric Gingival Crevicular Fluid FTIR
Measurements: Partition Stability and Sensitivity to Sex-Associated and
Clean-Strip-Aligned Variation

This directory is the R project root. The supplied final correction is applied
in full. The sole code change from the previously checked snapshot is:
  R/07_tables_figures.R: variable = colnames(X_snv),
in the Table S7 export. All original scientific data/model outputs are retained.

The author reports that R07 and R08 were rerun and both corresponding validators
passed. The delivered S7 CSV has 793 rows and five columns in this order:
  wavenumber, variable, PC1, PC2, PC3

Current update/rebuild route, only if needed later (not executed by this audit):
  source("R/07_tables_figures.R")
  source("tests/validate_tables_figures_outputs.R")
  source("R/08_render_reproducibility_report.R")
  source("tests/validate_reproducibility_report.R")

Do not rerun scripts 01–06 for this correction. They and all frozen outputs
are included unchanged to preserve the analytical provenance. Script 01 needs
the original raw OPUS.zip, which is not included in this processed-input package.
The four original unused-strip DPT controls and the processed source spectra
are retained. No end-to-end raw-input reproduction is claimed.

For the R/Quarto environment and package versions, see the existing sessionInfo
files in tables/ and data/pediatric/processed/. R07 checks its own required
packages; R08 additionally requires a Quarto executable on PATH. No new
environment lockfile, package installation, or analytical branch was invented.

Published figures: Figure1–Figure6 and FigureS1–FigureS4 (PNG and PDF).
The other 12 PNG files are active diagnostic outputs of the frozen analytical
modules, not retired publication figures. Keep them for provenance.
Rendered self-contained report: analysis/revision_report.html
Source: analysis/revision_report.qmd
Main tables: tables/Table1_..., Table2_..., Table3_...
Supplementary CSVs: supplement/TableS1_... through TableS7_...

Original output manifests: tables/07_output_manifest.csv (30 entries) and
tables/08_report_output_manifest.csv (5 entries). All paths, byte sizes and MD5
hashes were checked. The root SHA-256 manifest describes this exact snapshot.
If files are intentionally regenerated later, the SHA-256 manifest must also
be regenerated for that new snapshot; PDF creation dates may alter file hashes
without changing rendered content.

No GitHub publication or journal submission has been performed.
