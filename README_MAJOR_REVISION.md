# Pediatric GCF FTIR — Major Revision reproducibility add-on

**Manuscript:** *Exploratory Chemometric Analysis of Pediatric Gingival Crevicular Fluid FTIR Measurements: Partition Stability and Sensitivity to Sex-Associated and Clean-Substrate-Aligned Variation*  
**Manuscript ID:** `applsci-4662734`  
**Version:** `v1.1.0-pediatric-major-revision`  
**Prepared:** 8 November 2026  

## Scope and provenance

This add-on provides the **new R analyses and machine-readable summaries developed for Major Revision**. It supplements the previously released, frozen baseline:

- Repository: https://github.com/Litvinova1984/applsci-gcf-ftir
- Baseline release: https://github.com/Litvinova1984/applsci-gcf-ftir/releases/tag/v1.0.0-pediatric-exploratory (15 October 2026)

The baseline remains the reference for the original n = 25, 793-variable matrix, PCA scores, Ward.D2 k = 3 partition (15/4/6), figures, and resampling results. This archive does not contain raw OPUS or clinical records, and it does **not** claim raw-to-processed reproducibility.

## Contents

| Path | Purpose |
| --- | --- |
| `R/09a_major_revision_preliminary_reproduction.R` | Age/sex summaries, age-adjusted PC1, frozen-score k = 2 comparison |
| `R/09_major_revision_sensitivity.R` | 4,999 sex-label permutations; quantified removed variation/overlap; alternative k and spectral-window analyses with initial separate RNG seeds |
| `R/09b_resampling_seed_audit.R` | **Author-verified seed-aligned final** resampling analyses and audit against the frozen original |
| `tables/D_alternative_windows_and_k_SEED_ALIGNED.csv` | Full **final** alternative-specification results; source of truth for S10 |
| `tables/TableS8_removed_variation_overlap.csv` | Supplementary Table S8 summary |
| `tables/TableS9_age_puberty_by_sex.csv` | Supplementary Table S9 summary |
| `tables/TableS10_analytical_specification_diagnostic_seed_aligned.csv` | Supplementary Table S10 (aligned, rounded) |
| `tables/TableS11_sex_label_permutations.csv` | Supplementary Table S11 summary |
| `RELEASE_NOTES_GITHUB.md` | Ready-to-copy GitHub release description |
| `MANIFEST_SHA256.csv` | Integrity fingerprints of packaged files, except the manifest itself |

## Reproduce from the baseline repository

1. Download or clone the **baseline repository**, including the processed `data/` and `R/00_utils.R` files. The add-on is not standalone.
2. Copy the three new scripts into the repository's `R/` directory (alongside scripts 00–08).
3. From the **repository root** in R/RStudio, execute:

```r
source("R/09a_major_revision_preliminary_reproduction.R")
source("R/09_major_revision_sensitivity.R")
source("R/09b_resampling_seed_audit.R")
```

The scripts write to `major_revision_outputs/`, `major_revision_outputs/preliminary/`, and `major_revision_outputs/seed_aligned/` respectively. Previously frozen files are not overwritten.

**CRITICAL:** `09_major_revision_sensitivity.R` intentionally uses different RNG seeds (`SEED + 1900`, `SEED + 1901`) for its preliminary alternative-specification diagnostics. **Table S10 and the final values must be taken from `09b_resampling_seed_audit.R`, i.e. `major_revision_outputs/seed_aligned/D_alternative_windows_and_k.csv`, not from the unaligned `major_revision_outputs/D_alternative_windows_and_k.csv`.** The two scripts are retained separately for full traceability.

The original baseline bootstrap/subsampling outputs remain frozen. The aligned final variant calculations use the original `SEED = 4489628`, resetting it before each subsampling and bootstrap draw pool, like `R/03_cluster_stability.R`.

## Author-reported validation and key checks

The author reported a successful local execution on the original processed dataset:

```text
Original subsampling mean ARI: 0.6418108439
Aligned subsampling mean ARI:  0.6418108439
Original bootstrap mean ARI:   0.5794283128
Aligned bootstrap mean ARI:    0.5794283128
PASS: both means match the original report to 1e-6.
```

Other reported results: 4,999 sex-label permutations, lower-tail p = 0.0784; removed centered spectral sum-of-squares 15.74568% (sex residualization) and 6.967528% (clean-reference direction); k = 2 vs frozen k = 3 ARI = 0.614792; fingerprint-only recalculated-SNV k = 3 ARI = 0.773816. Detailed final metrics are in `tables/`.



## Interpretation boundaries

- The k = 3 solution is a descriptive reference, **not** a biological phenotype classification.
- The removed-variation percentages are properties of mathematical transformations, **not** isolated biological or substrate fractions.
- Variable names such as `clean_strip` in historical code refer to measurements of **Meta Biomed endodontic absorbent paper points (ISO 20)**, not flat periodontal strips.
- The original raw OPUS files, source-linked clinical records, and individual technical acquisitions are not distributed.

## Integrity

`MANIFEST_SHA256.csv` lists the SHA-256 digest of each packaged file except the manifest. 
