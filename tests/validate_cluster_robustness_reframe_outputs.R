# ============================================================
# validate_cluster_robustness_reframe_outputs.R
# Numerical and scope guardrails for R/03_cluster_stability.R
# ============================================================
# Run after:
#   source("R/03_cluster_stability.R")
# from the repository root.
# ============================================================

required_files <- c(
  file.path("tables", "cluster_PC_count_sensitivity.csv"),
  file.path("tables", "cluster_number_sensitivity_k2_k6.csv"),
  file.path("tables", "baseline_silhouette_by_participant.csv"),
  file.path("tables", "baseline_silhouette_by_cluster.csv"),
  file.path("tables", "cluster_LOO_summary.csv"),
  file.path("tables", "cluster_subsampling_summary.csv"),
  file.path("tables", "cluster_consensus_summary.csv"),
  file.path("tables", "cluster_bootstrap_summary.csv"),
  file.path("tables", "cluster_bootstrap_Jaccard_summary.csv"),
  file.path("tables", "cluster_robustness_summary.csv"),
  file.path(
    "data",
    "pediatric",
    "processed",
    "cluster_robustness_objects.rds"
  )
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing output file(s). Run R/03_cluster_stability.R first:\n",
    paste(missing_files, collapse = "\n")
  )
}

legacy_files <- c(
  file.path("tables", "baseline_silhouette_full_snv_by_participant.csv"),
  file.path("tables", "baseline_internal_characterization_full_snv.csv"),
  file.path(
    "data",
    "pediatric",
    "processed",
    "baseline_internal_characterization_full_snv.rds"
  )
)

if (any(file.exists(legacy_files))) {
  stop(
    "Legacy full-SNV/PERMANOVA output(s) remain after the reframe:\n",
    paste(legacy_files[file.exists(legacy_files)], collapse = "\n")
  )
}

pc <- read.csv(
  file.path("tables", "cluster_PC_count_sensitivity.csv"),
  check.names = FALSE
)

k_tbl <- read.csv(
  file.path("tables", "cluster_number_sensitivity_k2_k6.csv"),
  check.names = FALSE
)

sil <- read.csv(
  file.path("tables", "baseline_silhouette_by_participant.csv"),
  check.names = FALSE
)

loo <- read.csv(
  file.path("tables", "cluster_LOO_summary.csv"),
  check.names = FALSE
)

subsample <- read.csv(
  file.path("tables", "cluster_subsampling_summary.csv"),
  check.names = FALSE
)

consensus <- read.csv(
  file.path("tables", "cluster_consensus_summary.csv"),
  check.names = FALSE
)

bootstrap <- read.csv(
  file.path("tables", "cluster_bootstrap_summary.csv"),
  check.names = FALSE
)

jaccard <- read.csv(
  file.path("tables", "cluster_bootstrap_Jaccard_summary.csv"),
  check.names = FALSE
)

summary_tbl <- read.csv(
  file.path("tables", "cluster_robustness_summary.csv"),
  check.names = FALSE
)

stopifnot(
  nrow(pc) == 6L,
  identical(as.integer(pc$n_PC), 2:7),
  nrow(k_tbl) == 5L,
  identical(as.integer(k_tbl$k), 2:6),
  nrow(sil) == 25L,
  sum(sil$silhouette_width < 0) == 0L,
  abs(mean(sil$silhouette_width) - 0.434) < 0.02,
  abs(loo$mean_ARI - 0.811) < 0.03,
  abs(loo$median_ARI - 1.000) < 1e-12,
  abs(loo$proportion_ARI_1 - 0.68) < 0.05,
  abs(subsample$mean_ARI - 0.647) < 0.04,
  abs(bootstrap$mean_ARI - 0.601) < 0.04,
  nrow(jaccard) == 3L
)

exact_pc <- as.integer(
  pc$n_PC[pc$identical_to_baseline]
)
stopifnot(identical(exact_pc, c(2L, 3L, 5L, 6L)))

within_mean <- consensus$mean_consensus[
  consensus$baseline_same_cluster %in% c(TRUE, "TRUE")
]
between_mean <- consensus$mean_consensus[
  consensus$baseline_same_cluster %in% c(FALSE, "FALSE")
]

stopifnot(
  length(within_mean) == 1L,
  length(between_mean) == 1L,
  within_mean > 0.70,
  between_mean < 0.20
)

expected_jaccard <- c(0.768, 0.547, 0.827)
observed_jaccard <- jaccard$mean_Jaccard[
  order(jaccard$baseline_cluster)
]
stopifnot(max(abs(observed_jaccard - expected_jaccard)) < 0.04)

forbidden_pattern <- "full.?SNV|PERMANOVA|betadisper"
summary_text <- paste(summary_tbl$section, summary_tbl$metric)
if (any(grepl(forbidden_pattern, summary_text, ignore.case = TRUE))) {
  stop("cluster_robustness_summary.csv still contains retired full-SNV inference.")
}

obj <- readRDS(
  file.path(
    "data",
    "pediatric",
    "processed",
    "cluster_robustness_objects.rds"
  )
)

forbidden_object_names <- c(
  "full_snv_silhouette_by_participant",
  "baseline_internal_characterization",
  "permanova_full_snv",
  "betadisper_full_snv",
  "betadisper_perm"
)

if (any(forbidden_object_names %in% names(obj))) {
  stop(
    "Retired full-SNV inference objects remain in cluster_robustness_objects.rds."
  )
}

script_text <- paste(
  readLines(file.path("R", "03_cluster_stability.R"), warn = FALSE),
  collapse = "\n"
)

# References inside the explicit stale-output cleanup are allowed; active
# analytical calls and dependencies are not.
if (grepl("vegan::|adonis2\\(|betadisper\\(|permutest\\(", script_text)) {
  stop("R/03_cluster_stability.R still contains retired inferential calls.")
}

cat("Cluster-robustness reframe validation passed.\n")
