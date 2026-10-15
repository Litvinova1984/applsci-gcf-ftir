# Applied Sciences exploratory resubmission — release checklist

**Scope:** final processed-input reproducibility snapshot, 15 October 2026. No new scientific analyses or changes to approved results.

## Before commit

- [ ] Current Git branch is `resubmission-2026-10-15`; old `main` is preserved until review/merge.
- [ ] `README.md` documents the exploratory paper, R02–08 processed-input rebuild, and raw OPUS restriction.
- [ ] `README.md` and this checklist contain no references to executable `R/run_all_public.R`, adult PCA or old D1–D5 analyses as active steps.
- [ ] Git-tracked files correspond to the checked `reproducibility/` snapshot, plus explicitly retained project setup and updated documentation.
- [ ] No `data/pediatric/raw/OPUS.zip`, private medical records, direct personal identifiers, personal local paths, or working ZIP archives are staged.
- [ ] `supplement/TableS7_PCA_eigenvector_loadings_PC1_PC3.csv` has 793 rows and five columns: `wavenumber,variable,PC1,PC2,PC3`.
- [ ] `R/07_tables_figures.R` exports column `variable`.
- [ ] Main Figure1–Figure6 and supplementary FigureS1–FigureS4 are present, in PNG and PDF.
- [ ] Seven supplementary CSVs, `analysis/revision_report.html`, tests, manifests and session information are present.
- [ ] Original archived `MANIFEST_SHA256.csv` still verifies its listed snapshot files; repository README/checklist need not be listed there.
- [ ] `git status` has been reviewed before staging and pushing.

## QA provenance

- Final package QA: **11 PASS, 0 FAIL, 2 NOT VERIFIED**.
- Author successfully reran R07 and R08 and the corresponding validators.
- Independent R execution by the assembler: **NOT VERIFIED**.
- Full raw OPUS-to-output reproduction: **NOT VERIFIED**; restricted raw OPUS is not public.

## Publication sequence

1. Stage and review changes; commit to the dedicated branch.
2. Push the branch to GitHub and review the complete file diff before merging to `main`.
3. Merge the approved exploratory version to `main` (if that is the intended release branch).
4. Create a **versioned GitHub Release** from the approved commit, with a descriptive tag and accurate scope notes. Do not claim full raw-input reproduction.
5. Obtain the exact release URL; add it consistently to the manuscript Data Availability Statement and Cover Letter.
6. Verify final publication documents and submit through the journal system only on explicit author approval.
