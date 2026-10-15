# ============================================================
# validate_sex_reframe_outputs.R
# Numerical guardrails for R/04_sex_sensitivity.R
# ============================================================

required_files <- c(
  file.path("tables", "sex_baseline_cluster_association.csv"),
  file.path("tables", "PC1_PC3_sex_associations_BH.csv"),
  file.path("tables", "sex_residualized_cluster_assignments.csv"),
  file.path("tables", "sex_residualized_clustering_summary.csv"),
  file.path("tables", "sex_sensitivity_summary.csv"),
  file.path(
    "data",
    "pediatric",
    "processed",
    "sex_sensitivity_objects.rds"
  )
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing output file(s). Run R/04_sex_sensitivity.R first:\n",
    paste(missing_files, collapse = "\n")
  )
}

sex_baseline <- read.csv(
  file.path("tables", "sex_baseline_cluster_association.csv"),
  check.names = FALSE
)

pc_results <- read.csv(
  file.path("tables", "PC1_PC3_sex_associations_BH.csv"),
  check.names = FALSE
)

assignments <- read.csv(
  file.path("tables", "sex_residualized_cluster_assignments.csv"),
  check.names = FALSE
)

summary_tbl <- read.csv(
  file.path("tables", "sex_residualized_clustering_summary.csv"),
  check.names = FALSE
)

stopifnot(
  nrow(sex_baseline) == 1L,
  nrow(pc_results) == 3L,
  identical(as.character(pc_results$PC), c("PC1", "PC2", "PC3")),
  nrow(assignments) == 25L,
  length(unique(assignments$participant_id)) == 25L,
  nrow(summary_tbl) == 1L
)

bh_expected <- stats::p.adjust(pc_results$p_value_raw, method = "BH")
stopifnot(max(abs(pc_results$p_value_BH - bh_expected)) < 1e-12)

# Frozen-result guardrails from the submitted analysis. The intervals
# allow for negligible numerical differences between R/package builds.
stopifnot(
  sex_baseline$fisher_p_value < 0.01,
  sex_baseline$cramers_v > 0.50,
  sex_baseline$cramers_v < 0.80,
  pc_results$p_value_raw[pc_results$PC == "PC1"] < 0.05,
  pc_results$p_value_BH[pc_results$PC == "PC1"] < 0.05,
  all(pc_results$p_value_BH[pc_results$PC %in% c("PC2", "PC3")] > 0.05),
  abs(summary_tbl$ARI_vs_baseline - 0.134) < 0.02,
  abs(summary_tbl$residualized_average_silhouette - 0.338) < 0.03,
  abs(
    summary_tbl$residualized_3PC_cumulative_variance_percent - 87.5
  ) < 0.5,
  summary_tbl$sex_cluster_fisher_p_after_residualization > 0.05,
  summary_tbl$max_abs_design_residual_crossproduct < 1e-8
)

cluster_sizes <- sort(
  as.integer(
    strsplit(summary_tbl$residualized_cluster_sizes, "/", fixed = TRUE)[[1]]
  )
)
stopifnot(identical(cluster_sizes, c(3L, 8L, 14L)))

object_bundle <- readRDS(
  file.path(
    "data",
    "pediatric",
    "processed",
    "sex_sensitivity_objects.rds"
  )
)

stopifnot(
  identical(
    object_bundle$settings$pc_multiplicity_adjustment,
    "Benjamini-Hochberg across PC1-PC3"
  )
)

cat("Sex-reframe validation passed.\n")
