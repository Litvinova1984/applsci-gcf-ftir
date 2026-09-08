# Run the public reproducibility workflow from the repository root.

scripts <- c(
  "R/01_validate_public_inputs.R",
  "R/02_pca_clustering.R",
  "R/03_cluster_stability.R",
  "R/04_sex_sensitivity.R",
  "R/05_strip_reproducibility.R",
  "R/06_spectral_domains.R",
  "R/07_adult_reference_public.R",
  "R/08_tables_figures.R"
)

for (script in scripts) {
  cat("\n\n>>> Running ", script, "\n", sep = "")
  source(script, echo = FALSE)
}

cat("\nPublic analytical workflow completed successfully.\n")
cat("Render the reproducibility report with:\n")
cat('  source("R/render_revision_report.R")\n')
