# ============================================================
# 08_render_reproducibility_report.R
#
# Render the pediatric-only Quarto reproducibility report.
# Run from the project root after scripts 01-07 and their validators.
# ============================================================

qmd <- file.path("analysis", "revision_report.qmd")
html <- file.path("analysis", "revision_report.html")
tables_dir <- "tables"

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(qmd)) {
  stop("Missing ", qmd, ".")
}

required_before_render <- c(
  file.path("tables", "07_output_manifest.csv"),
  file.path("tables", "Table1_pediatric_cohort_characteristics.csv"),
  file.path("tables", "Table2_exploratory_partition_stability.csv"),
  file.path("tables", "Table3_sex_and_substrate_sensitivity.csv"),
  file.path("figures", "Figure1_analytical_workflow.png"),
  file.path("figures", "Figure2_spectral_structure_and_PCA_profiles.png"),
  file.path("figures", "Figure3_partition_selection_and_stability.png"),
  file.path("figures", "Figure4_sex_associated_variation.png"),
  file.path("figures", "Figure5_clean_strip_reproducibility_and_sensitivity.png"),
  file.path("figures", "Figure6_partition_assignment_transitions.png")
)

missing_before_render <- required_before_render[!file.exists(required_before_render)]
if (length(missing_before_render) > 0L) {
  stop(
    "Required analytical outputs are missing. Run scripts 01-07 first:\n",
    paste(missing_before_render, collapse = "\n")
  )
}

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  stop(
    "Quarto executable was not found on PATH. ",
    "Open analysis/revision_report.qmd in RStudio/Posit and click Render, ",
    "or install Quarto and rerun this script."
  )
}

if (file.exists(html)) unlink(html, force = TRUE)

status <- system2(
  quarto_bin,
  c("render", qmd, "--to", "html")
)

if (!identical(status, 0L)) {
  stop("Quarto render failed with exit status ", status, ".")
}

required_after_render <- c(
  qmd,
  html,
  file.path(tables_dir, "08_reproducibility_report_key_numbers.csv"),
  file.path(tables_dir, "08_reproducibility_report_scope_checks.csv")
)

missing_after_render <- required_after_render[!file.exists(required_after_render)]
if (length(missing_after_render) > 0L) {
  stop(
    "Report rendering did not create all expected outputs:\n",
    paste(missing_after_render, collapse = "\n")
  )
}

if (file.info(html)$size <= 0) {
  stop("Rendered HTML is empty.")
}

quarto_version <- tryCatch(
  system2(quarto_bin, "--version", stdout = TRUE, stderr = TRUE),
  error = function(e) paste("Unable to query Quarto version:", conditionMessage(e))
)

session_path <- file.path(tables_dir, "sessionInfo_reproducibility_report.txt")
capture.output(
  {
    cat("Quarto executable: ", quarto_bin, "\n", sep = "")
    cat("Quarto version: ", paste(quarto_version, collapse = " "), "\n\n", sep = "")
    print(sessionInfo())
  },
  file = session_path
)

report_paths <- c(required_after_render, session_path)
report_manifest <- data.frame(
  category = c("source", "rendered_report", "audit_table", "audit_table", "session_info"),
  path = report_paths,
  exists = file.exists(report_paths),
  size_bytes = as.numeric(file.info(report_paths)$size),
  md5 = unname(tools::md5sum(report_paths)),
  stringsAsFactors = FALSE
)

if (!all(report_manifest$exists) || any(report_manifest$size_bytes <= 0)) {
  stop("One or more report outputs are missing or empty.")
}

utils::write.csv(
  report_manifest,
  file.path(tables_dir, "08_report_output_manifest.csv"),
  row.names = FALSE,
  na = ""
)

cat("\n========================================\n")
cat("PEDIATRIC REPRODUCIBILITY REPORT RENDERED\n")
cat("========================================\n")
cat("HTML: ", normalizePath(html, mustWork = TRUE), "\n", sep = "")
cat(
  "Key-number snapshot: ",
  normalizePath(
    file.path(tables_dir, "08_reproducibility_report_key_numbers.csv"),
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat("Manifest: tables/08_report_output_manifest.csv\n")
cat("========================================\n")
