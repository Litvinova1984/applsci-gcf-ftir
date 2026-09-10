# Candidate FTIR-Derived Spectral Phenotypes of Pediatric Gingival Crevicular Fluid Identified by Chemometric Analysis

Public reproducibility repository for the revised *Applied Sciences* manuscript.

## Scope

This repository reproduces the revised analyses from de-identified processed pediatric FTIR spectra, de-identified participant-level pediatric metadata, four independent clean-strip technical-reference spectra, a de-identified processed adult reference spectral matrix, and frozen **aggregate** adult clinical/axis summaries.

The public workflow does **not** distribute source-linked clinical records, pediatric raw OPUS archives, adult source labels, individual-level adult clinical metadata, or private audit files. The adult reference is used only for exploratory domain-level spectroscopic comparison; pediatric and adult spectra are never pooled into a common PCA model and the adult data are not used for pediatric classification or individual prediction.

## Frozen analytical specification

- Pediatric active cohort: **n = 25**
- Retained spectral variables: **793**
- Analytical windows: **3400–2800 and 1800–870 cm⁻¹**
- Preprocessing: row-wise standard normal variate (SNV); no additional feature scaling
- Baseline PCA: first **3 PCs**, cumulative variance **89.282%**
- Baseline clustering: Euclidean distance + Ward.D2 hierarchical clustering
- Baseline cluster sizes: **15 / 4 / 6**
- Fixed random seed: **4489628**
- Adult reference: **18 × 793** harmonized SNV matrix
- Historical-to-harmonized adult axis correspondence: **PC1→PC2; PC2→PC1; PC4→PC4; PC5→PC5**
- Adult domain map: **D1/D2→PC2; D3/D4→PC4; D5→PC5**

## Repository structure

```text
.
├── R/
│   ├── 00_utils.R
│   ├── 01_validate_public_inputs.R
│   ├── 02_pca_clustering.R
│   ├── 03_cluster_stability.R
│   ├── 04_sex_sensitivity.R
│   ├── 05_strip_reproducibility.R
│   ├── 06_spectral_domains.R
│   ├── 07_adult_reference_public.R
│   ├── 08_tables_figures.R
│   ├── render_revision_report.R
│   └── run_all_public.R
├── analysis/
│   └── revision_report.qmd
├── data/
│   ├── pediatric/
│   │   ├── pediatric_metadata.csv
│   │   └── processed/
│   ├── adult/
│   │   ├── processed/
│   │   └── aggregate/
│   └── technical_controls/
├── figures/
│   └── source_assets/
├── tables/
├── supplement/
└── renv.lock
```

## Public inputs

### Pediatric

```text
data/pediatric/pediatric_metadata.csv
data/pediatric/processed/pediatric_X_active_absorbance.rds
data/pediatric/processed/pediatric_X_snv.rds
data/pediatric/processed/pediatric_wavenumbers.rds
data/pediatric/processed/pediatric_source_spectra_long.rds
```

The public long-format spectral object contains the analytical spectra required for the supplementary pediatric references and archived clean-strip reference. Its source labels are study-local generic analytical filenames only (e.g., numeric `.spc` labels and `etalon/G/N-*` labels); no patient names, dates, filesystem paths, or clinical-record identifiers are included.

### Technical clean-strip references

```text
data/technical_controls/etalon_01.dpt
data/technical_controls/etalon_02.dpt
data/technical_controls/etalon_03.dpt
data/technical_controls/etalon_04.dpt
```

These are four independent unused strips from the same production batch.

### Adult reference

```text
data/adult/processed/adult_reference_X_snv.rds
data/adult/processed/adult_reference_wavenumbers.rds
data/adult/aggregate/*.csv
```

The `aggregate/` directory contains only frozen group-level/axis-level outputs needed to document the historical-axis matching and adult clinical context. Individual-level adult clinical metadata and adult source labels are not included.

`R/07_adult_reference_public.R` independently recomputes the harmonized adult PCA and D1–D5 loading summaries from the public 18 × 793 spectral matrix, verifies them against frozen aggregate loading/domain summaries, and then creates a minimal aggregate-only R object consumed by script 08 and the Quarto report.

### Curated workflow figure

```text
figures/source_assets/Figure1_analytical_workflow.png
```

Figure 1 is a curated workflow schematic and is copied into the publication-figure directory by script 08.

## Running the analysis

Run from the repository root in R/RStudio:

```r
install.packages("renv")
renv::restore()
source("R/run_all_public.R")
```

Script 08 regenerates the final tables and figures, including the final color Supplementary Figures S1–S5. It also writes `tables/08_output_manifest.csv`.

Then render the reproducibility report:

```r
source("R/render_revision_report.R")
```

Quarto must be installed separately and available on the system path. The rendered report is written to `analysis/revision_report.html`.

## Main output set

- Table 1
- Figures 1–7
- Supplementary Tables S1–S5/S5b
- Supplementary Figures S1–S5
- aggregate adult reference/loading exports
- `tables/08_output_manifest.csv`
- `analysis/revision_report.html`

## Interpretation and validation scope

The three pediatric groups are treated as **candidate FTIR-derived spectral phenotypes**, not validated clinical classes or pathophysiological endotypes. PLS-DA/VIP is used only for post-clustering spectral characterization, not as independent validation of the cluster solution.

Sex residualization, clean-strip residualization, and D1-region exclusion are sensitivity analyses, not alternative corrected ground-truth partitions.

The archived adult dataset is an exploratory spectroscopic reference only and is not a diagnostic or predictive model for pediatric GCF.

## Data privacy

The repository intentionally excludes:

- raw pediatric OPUS archives;
- source-linked pediatric clinical records;
- potentially identifying provenance information;
- adult source labels;
- individual-level adult clinical metadata;
- private audit files.

The pediatric participant identifiers in the public analytical data are study-local de-identified IDs.

## Environment

The R package environment is recorded in `renv.lock`. Session information produced by the analysis modules is written to `tables/`.
