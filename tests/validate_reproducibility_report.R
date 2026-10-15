# ============================================================
# validate_reproducibility_report.R
#
# Validate the rendered pediatric-only Quarto report.
# Run from the project root after R/08_render_reproducibility_report.R.
# ============================================================

required_files <- c(
  file.path("analysis", "revision_report.qmd"),
  file.path("analysis", "revision_report.html"),
  file.path("tables", "08_reproducibility_report_key_numbers.csv"),
  file.path("tables", "08_reproducibility_report_scope_checks.csv"),
  file.path("tables", "08_report_output_manifest.csv"),
  file.path("tables", "sessionInfo_reproducibility_report.txt")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required report output(s) are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

sizes <- file.info(required_files)$size
if (any(!is.finite(sizes)) || any(sizes <= 0)) {
  stop("At least one required report output is empty or unreadable.")
}

key_numbers <- utils::read.csv(
  file.path("tables", "08_reproducibility_report_key_numbers.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

scope_checks <- utils::read.csv(
  file.path("tables", "08_reproducibility_report_scope_checks.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manifest <- utils::read.csv(
  file.path("tables", "08_report_output_manifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(key_numbers) != 39L) {
  stop("Expected 39 key-number rows; observed ", nrow(key_numbers), ".")
}

if (nrow(scope_checks) != 6L || !all(as.logical(scope_checks$passed))) {
  stop("One or more report scope checks failed.")
}

if (nrow(manifest) != 5L || !all(as.logical(manifest$exists))) {
  stop("Unexpected report manifest structure.")
}

get_value <- function(section, metric) {
  z <- key_numbers$value_numeric[
    key_numbers$section == section & key_numbers$metric == metric
  ]
  if (length(z) != 1L) {
    stop("Expected one key-number row for ", section, " / ", metric, ".")
  }
  as.numeric(z[[1]])
}

if (get_value("Data", "Participants") != 25) {
  stop("Participant count changed unexpectedly.")
}
if (get_value("Data", "Spectral variables") != 793) {
  stop("Spectral-variable count changed unexpectedly.")
}
if (abs(get_value("Baseline PCA", "Cumulative variance PC1-PC3 (%)") - 89.282) > 0.05) {
  stop("Unexpected cumulative variance for PC1-PC3.")
}
if (!identical(
  as.integer(c(
    get_value("Baseline partition", "Group 1 n"),
    get_value("Baseline partition", "Group 2 n"),
    get_value("Baseline partition", "Group 3 n")
  )),
  c(15L, 4L, 6L)
)) {
  stop("Baseline group sizes changed unexpectedly.")
}
if (abs(get_value("Sex-associated variation", "Sex-residualized ARI vs baseline") - 0.134) > 0.02) {
  stop("Unexpected sex-residualized ARI.")
}
if (abs(get_value("Clean-strip analysis", "Residualized ARI vs baseline") - 0.407) > 0.02) {
  stop("Unexpected clean-strip-residualized ARI.")
}
if (get_value("Clean-strip analysis", "Minimum pairwise Pearson r") < 0.9999) {
  stop("Unexpectedly low clean-strip reproducibility.")
}

qmd_text <- paste(readLines(file.path("analysis", "revision_report.qmd"), warn = FALSE), collapse = "\n")
html_text <- paste(readLines(file.path("analysis", "revision_report.html"), warn = FALSE), collapse = "\n")

required_headings <- c(
  "Scope and computational contract",
  "Spectral structure and baseline exploratory partition",
  "Partition selection and internal stability",
  "Sex-associated spectral variation",
  "Clean-strip technical reproducibility and sensitivity",
  "Participant-level partition transitions",
  "Output manifest and machine-readable audit"
)

if (any(!vapply(required_headings, grepl, logical(1), x = html_text, fixed = TRUE))) {
  stop("The rendered report is missing one or more required sections.")
}

required_figure_names <- c(
  "Figure1_analytical_workflow.png",
  "Figure2_spectral_structure_and_PCA_profiles.png",
  "Figure3_partition_selection_and_stability.png",
  "Figure4_sex_associated_variation.png",
  "Figure5_clean_strip_reproducibility_and_sensitivity.png",
  "Figure6_partition_assignment_transitions.png"
)

if (any(!vapply(required_figure_names, grepl, logical(1), x = qmd_text, fixed = TRUE))) {
  stop("The report source is missing one or more main figures.")
}

retired_visible_phrases <- c(
  "candidate FTIR-derived spectral phenotype",
  "early or transitional disease-related state",
  "adult-reference mapping",
  "PLS-DA",
  "VIP-enriched"
)

if (any(vapply(retired_visible_phrases, grepl, logical(1), x = html_text, fixed = TRUE))) {
  stop("A retired interpretive phrase remains in the rendered report.")
}

if (file.info(file.path("analysis", "revision_report.html"))$size < 10000) {
  stop("Rendered report is unexpectedly small.")
}

cat("\n========================================\n")
cat("REPRODUCIBILITY REPORT VALIDATION PASSED\n")
cat("========================================\n")
cat("Key-number rows: ", nrow(key_numbers), "\n", sep = "")
cat("Scope checks: ", nrow(scope_checks), "\n", sep = "")
cat("Report manifest rows: ", nrow(manifest), "\n", sep = "")
cat("HTML size: ", file.info(file.path("analysis", "revision_report.html"))$size, " bytes\n", sep = "")
cat("========================================\n")
