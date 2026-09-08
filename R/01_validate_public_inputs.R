# ============================================================
# 01_validate_public_inputs.R
#
# Pediatric GCF FTIR — validation of de-identified public inputs
#
# The public repository starts from processed, de-identified spectra.
# Raw OPUS archives and potentially identifying provenance are not distributed.
# The long-format spectral object may retain study-local generic analytical
# filenames (numeric .spc labels and etalon/G/N-* reference labels only).
#
# Required public inputs:
#   data/pediatric/processed/pediatric_X_active_absorbance.rds
#   data/pediatric/processed/pediatric_X_snv.rds
#   data/pediatric/processed/pediatric_wavenumbers.rds
#   data/pediatric/processed/pediatric_source_spectra_long.rds
#   data/pediatric/pediatric_metadata.csv
#
# This script performs integrity checks only; it introduces no new
# preprocessing or inferential analysis.
# ============================================================

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}
source(file.path("R", "00_utils.R"))

required_pkgs <- c("readr", "dplyr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

processed_dir <- file.path("data", "pediatric", "processed")
metadata_file <- file.path("data", "pediatric", "pediatric_metadata.csv")

required_files <- c(
  file.path(processed_dir, "pediatric_X_active_absorbance.rds"),
  file.path(processed_dir, "pediatric_X_snv.rds"),
  file.path(processed_dir, "pediatric_wavenumbers.rds"),
  file.path(processed_dir, "pediatric_source_spectra_long.rds"),
  metadata_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing public input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

X_raw <- as.matrix(readRDS(file.path(processed_dir, "pediatric_X_active_absorbance.rds")))
X_snv <- as.matrix(readRDS(file.path(processed_dir, "pediatric_X_snv.rds")))
wn <- as.numeric(readRDS(file.path(processed_dir, "pediatric_wavenumbers.rds")))
spectra_long <- readRDS(file.path(processed_dir, "pediatric_source_spectra_long.rds"))
metadata <- readr::read_csv(metadata_file, show_col_types = FALSE)

if (!identical(dim(X_raw), c(25L, 793L))) stop("Expected X_raw = 25 x 793.")
if (!identical(dim(X_snv), c(25L, 793L))) stop("Expected X_snv = 25 x 793.")
if (length(wn) != 793L) stop("Expected 793 retained wavenumbers.")
if (!identical(dim(X_raw), dim(X_snv))) stop("Raw and SNV matrices differ in dimensions.")
if (any(!is.finite(X_raw)) || any(!is.finite(X_snv)) || any(!is.finite(wn))) {
  stop("Non-finite values detected in public spectral inputs.")
}

if (is.null(rownames(X_raw)) || is.null(rownames(X_snv))) {
  stop("Pediatric matrices must contain de-identified participant IDs as row names.")
}
if (!identical(rownames(X_raw), rownames(X_snv))) {
  stop("Raw and SNV matrices have different participant ordering.")
}
if (is.null(colnames(X_snv))) stop("X_snv must have wavenumber column names.")

wn_from_names <- extract_wavenumber(colnames(X_snv))
if (anyNA(wn_from_names)) stop("Could not parse wavenumbers from X_snv column names.")
if (max(abs(wn_from_names - wn)) > 0.01) {
  stop("Stored precise wavenumbers and matrix labels differ unexpectedly.")
}

if (max(abs(rowMeans(X_snv))) > 1e-10) {
  stop("SNV check failed: row means are not approximately zero.")
}
if (max(abs(apply(X_snv, 1, stats::sd) - 1)) > 1e-10) {
  stop("SNV check failed: row SDs are not approximately one.")
}

required_long_cols <- c("sample_id", "wavenumber", "absorbance")
if (!all(required_long_cols %in% names(spectra_long))) {
  stop(
    "pediatric_source_spectra_long.rds must contain: ",
    paste(required_long_cols, collapse = ", ")
  )
}
if (any(!is.finite(spectra_long$wavenumber)) ||
    any(!is.finite(spectra_long$absorbance))) {
  stop("Non-finite values detected in supplementary/reference spectra.")
}

expected_supp <- c("etalon", "G", "N-1", "N-2", "N")
if (!all(expected_supp %in% unique(as.character(spectra_long$sample_id)))) {
  stop("Expected supplementary/reference spectra are missing.")
}

required_meta <- c(
  "participant_id", "age", "sex", "pubertal_development",
  "endocrine_status", "uctd_status", "periodontal_status"
)
if (!all(required_meta %in% names(metadata))) {
  stop("Public pediatric metadata do not contain the expected columns.")
}
if (nrow(metadata) != 25L) stop("Expected 25 pediatric metadata rows.")

cat("\n========================================\n")
cat("PUBLIC PEDIATRIC INPUTS VALIDATED\n")
cat("========================================\n")
cat("Active pediatric spectra: 25\n")
cat("Retained spectral variables: 793\n")
cat("Analytical windows: 3400-2800 + 1800-870 cm^-1\n")
cat("SNV row-wise normalization checks: passed\n")
cat("Supplementary/reference spectra: present\n")
cat("De-identified participant metadata: 25 rows\n")
cat("========================================\n")
