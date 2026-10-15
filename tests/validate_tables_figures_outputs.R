# ============================================================
# validate_tables_figures_outputs.R
#
# Structural and scope validation for R/07_tables_figures.R
# Run from the project root after source("R/07_tables_figures.R").
# ============================================================

required_files <- c(
  file.path("tables", "Table1_pediatric_cohort_characteristics.csv"),
  file.path("tables", "Table2_exploratory_partition_stability.csv"),
  file.path("tables", "Table3_sex_and_substrate_sensitivity.csv"),
  file.path("tables", "07_output_manifest.csv"),
  file.path("supplement", "TableS1_cluster_stability_detailed.csv"),
  file.path("supplement", "TableS2_PC1_PC3_sex_associations_BH.csv"),
  file.path("supplement", "TableS3_clean_strip_technical_reproducibility.csv"),
  file.path("supplement", "TableS4_clean_strip_direction_sensitivity.csv"),
  file.path("supplement", "TableS5_participant_assignment_transitions.csv"),
  file.path("supplement", "TableS6_pairwise_partition_agreement.csv"),
  file.path("supplement", "TableS7_PCA_eigenvector_loadings_PC1_PC3.csv"),
  file.path("figures", "Figure1_analytical_workflow.png"),
  file.path("figures", "Figure1_analytical_workflow.pdf"),
  file.path("figures", "Figure2_spectral_structure_and_PCA_profiles.png"),
  file.path("figures", "Figure2_spectral_structure_and_PCA_profiles.pdf"),
  file.path("figures", "Figure3_partition_selection_and_stability.png"),
  file.path("figures", "Figure3_partition_selection_and_stability.pdf"),
  file.path("figures", "Figure4_sex_associated_variation.png"),
  file.path("figures", "Figure4_sex_associated_variation.pdf"),
  file.path("figures", "Figure5_clean_strip_reproducibility_and_sensitivity.png"),
  file.path("figures", "Figure5_clean_strip_reproducibility_and_sensitivity.pdf"),
  file.path("figures", "Figure6_partition_assignment_transitions.png"),
  file.path("figures", "Figure6_partition_assignment_transitions.pdf"),
  file.path("figures", "FigureS1_consensus_heatmap.png"),
  file.path("figures", "FigureS1_consensus_heatmap.pdf"),
  file.path("figures", "FigureS2_bootstrap_group_Jaccard.png"),
  file.path("figures", "FigureS2_bootstrap_group_Jaccard.pdf"),
  file.path("figures", "FigureS3_descriptive_group_spectra.png"),
  file.path("figures", "FigureS3_descriptive_group_spectra.pdf"),
  file.path("figures", "FigureS4_clean_strip_direction_invariance.png"),
  file.path("figures", "FigureS4_clean_strip_direction_invariance.pdf")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Required table/figure output(s) are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

file_sizes <- file.info(required_files)$size
if (any(!is.finite(file_sizes)) || any(file_sizes <= 0)) {
  stop("At least one required output file is empty or unreadable.")
}

manifest <- utils::read.csv(
  file.path("tables", "07_output_manifest.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_manifest_columns <- c(
  "category",
  "path",
  "exists",
  "size_bytes",
  "md5"
)

if (!all(required_manifest_columns %in% names(manifest))) {
  stop("The output manifest does not contain the expected columns.")
}

if (nrow(manifest) != 30L) {
  stop("Expected 30 manifest rows; observed ", nrow(manifest), ".")
}

if (!all(as.logical(manifest$exists))) {
  stop("The output manifest contains missing files.")
}

expected_category_counts <- c(
  main_table = 3L,
  supplement_table = 7L,
  main_figure = 12L,
  supplement_figure = 8L
)

observed_category_counts <- table(manifest$category)
for (nm in names(expected_category_counts)) {
  if (
    !(nm %in% names(observed_category_counts)) ||
      observed_category_counts[[nm]] != expected_category_counts[[nm]]
  ) {
    stop("Unexpected manifest count for category: ", nm)
  }
}

forbidden_manifest_pattern <- paste(
  c(
    "adult",
    "orientation",
    "domain",
    "phenotype",
    "caries",
    "PLS",
    "VIP",
    "PERMANOVA",
    "betadisper"
  ),
  collapse = "|"
)

if (any(grepl(forbidden_manifest_pattern, manifest$path, ignore.case = TRUE))) {
  stop("A retired or out-of-scope output appears in the new manifest.")
}

Table1 <- utils::read.csv(
  file.path("tables", "Table1_pediatric_cohort_characteristics.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
Table2 <- utils::read.csv(
  file.path("tables", "Table2_exploratory_partition_stability.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
Table3 <- utils::read.csv(
  file.path("tables", "Table3_sex_and_substrate_sensitivity.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(Table1) != 15L) {
  stop("Table 1 should contain 15 reporting rows.")
}

if (nrow(Table2) != 22L) {
  stop("Table 2 should contain 22 reporting rows.")
}

if (nrow(Table3) != 21L) {
  stop("Table 3 should contain 21 reporting rows.")
}

if (any(grepl("PERMANOVA|betadisper", Table2$Metric, ignore.case = TRUE))) {
  stop("Circular full-space inferential metrics remain in Table 2.")
}

if (!all(c("PC1 Wilcoxon BH-adjusted p", "ARI vs baseline") %in% Table3$Metric)) {
  stop("Table 3 does not contain the required sensitivity summaries.")
}

transition_table <- utils::read.csv(
  file.path("supplement", "TableS5_participant_assignment_transitions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(transition_table) != 25L) {
  stop("The participant transition table must contain 25 rows.")
}

if (anyDuplicated(transition_table$participant_id) > 0L) {
  stop("Duplicate participant IDs detected in the transition table.")
}

loading_table <- utils::read.csv(
  file.path("supplement", "TableS7_PCA_eigenvector_loadings_PC1_PC3.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(loading_table) != 793L) {
  stop("The PCA loading export must contain 793 spectral variables.")
}

if (!all(c("wavenumber", "PC1", "PC2", "PC3") %in% names(loading_table))) {
  stop("The PCA loading export does not contain the expected columns.")
}

retired_paths <- c(
  file.path("tables", "Table2_adult_informed_orientation_score.csv"),
  file.path("supplement", "TableS1_adult_directional_map.csv"),
  file.path("figures", "Figure5_D1_D5_domain_distributions.png"),
  file.path("figures", "Figure7_signed_domain_heatmap.png"),
  file.path("figures", "Figure8_adult_informed_orientation_score.png"),
  file.path("figures", "Figure9_control_proximity_vs_adult_orientation.png"),
  file.path("figures", "FigureS3_substrate_D1_sensitivity.png"),
  file.path("figures", "FigureS4_adult_PC1_PC4_PC5_loading_profiles.png")
)

if (any(file.exists(retired_paths))) {
  stop(
    "One or more retired publication outputs still exist:\n",
    paste(retired_paths[file.exists(retired_paths)], collapse = "\n")
  )
}

cat("\n========================================\n")
cat("TABLE / FIGURE OUTPUT VALIDATION PASSED\n")
cat("========================================\n")
cat("Required outputs checked: ", length(required_files), "\n", sep = "")
cat("Manifest rows checked: ", nrow(manifest), "\n", sep = "")
cat("Participant transition rows: ", nrow(transition_table), "\n", sep = "")
cat("PCA loading rows: ", nrow(loading_table), "\n", sep = "")
cat("========================================\n")
