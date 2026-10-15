# ============================================================
# 06_assignment_transitions.R
#
# Pediatric GCF FTIR exploratory reframe
# Consolidated participant-level partition transitions
#
# Run order:
#   01_import_preprocess.R
#   02_pca_clustering.R
#   03_cluster_stability.R
#   04_sex_sensitivity.R
#   05_strip_reproducibility.R
#   06_assignment_transitions.R
#
# PURPOSE
# -------
# This script combines the three participant partitions retained in
# the reframed Applied Sciences analysis:
#
#   - baseline PCA/Ward.D2 partition;
#   - sex-residualized partition;
#   - clean-strip-direction-residualized partition.
#
# The two sensitivity partitions have already been relabelled to
# maximize agreement with the baseline labels. The combined table is
# therefore intended for readable participant-level transition plots
# and counts. Adjusted Rand indices (ARI) remain label-invariant and
# are calculated independently from the displayed labels.
# ============================================================


# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "00_utils.R"))

required_pkgs <- c(
  "mclust",
  "tibble",
  "dplyr",
  "tidyr",
  "readr"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

processed_dir <- file.path("data", "pediatric", "processed")
tables_dir <- "tables"

baseline_cluster_file <- file.path(
  processed_dir,
  "pediatric_clusters_baseline.rds"
)

sex_assignment_file <- file.path(
  tables_dir,
  "sex_residualized_cluster_assignments.csv"
)

strip_assignment_file <- file.path(
  tables_dir,
  "strip_primary_assignment_transition.csv"
)

sex_summary_file <- file.path(
  tables_dir,
  "sex_residualized_clustering_summary.csv"
)

strip_summary_file <- file.path(
  tables_dir,
  "strip_primary_sensitivity_summary.csv"
)

required_files <- c(
  baseline_cluster_file,
  sex_assignment_file,
  strip_assignment_file,
  sex_summary_file,
  strip_summary_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required output file(s) are missing. Run scripts 02, 04, and 05 first:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------
# 2. Load and validate assignments
# ------------------------------------------------------------

baseline_raw <- readRDS(baseline_cluster_file)

if (is.null(names(baseline_raw))) {
  stop("The baseline cluster vector must have participant IDs as names.")
}

baseline <- if (is.factor(baseline_raw)) {
  as.integer(as.character(baseline_raw))
} else {
  as.integer(baseline_raw)
}
names(baseline) <- names(baseline_raw)

ids <- names(baseline)

if (length(ids) != 25L || anyDuplicated(ids) > 0L) {
  stop("Expected 25 unique participant IDs in the baseline partition.")
}

sex_assignments <- readr::read_csv(
  sex_assignment_file,
  show_col_types = FALSE
) |>
  dplyr::mutate(participant_id = as.character(participant_id))

strip_assignments <- readr::read_csv(
  strip_assignment_file,
  show_col_types = FALSE
) |>
  dplyr::mutate(participant_id = as.character(participant_id))

required_sex_columns <- c(
  "participant_id",
  "sex",
  "baseline_cluster",
  "sex_residualized_cluster",
  "assignment_changed"
)

required_strip_columns <- c(
  "participant_id",
  "baseline_cluster",
  "clean_strip_residualized_cluster",
  "assignment_changed"
)

if (!all(required_sex_columns %in% names(sex_assignments))) {
  stop("The sex-assignment file does not contain the expected columns.")
}

if (!all(required_strip_columns %in% names(strip_assignments))) {
  stop("The clean-strip assignment file does not contain the expected columns.")
}

if (
  anyDuplicated(sex_assignments$participant_id) > 0L ||
    anyDuplicated(strip_assignments$participant_id) > 0L
) {
  stop("Duplicate participant IDs detected in sensitivity assignments.")
}

if (
  !setequal(ids, sex_assignments$participant_id) ||
    !setequal(ids, strip_assignments$participant_id)
) {
  stop("Sensitivity-assignment participant IDs do not match baseline IDs.")
}

sex_assignments <- sex_assignments[
  match(ids, sex_assignments$participant_id),
  ,
  drop = FALSE
]

strip_assignments <- strip_assignments[
  match(ids, strip_assignments$participant_id),
  ,
  drop = FALSE
]

if (
  !identical(
    as.integer(sex_assignments$baseline_cluster),
    as.integer(baseline[ids])
  ) ||
    !identical(
      as.integer(strip_assignments$baseline_cluster),
      as.integer(baseline[ids])
    )
) {
  stop("Baseline labels disagree across the participant-level files.")
}


# ------------------------------------------------------------
# 3. Combined participant-level transition table
# ------------------------------------------------------------

assignment_transitions <- tibble::tibble(
  participant_id = ids,
  sex = as.character(sex_assignments$sex),
  baseline_cluster = as.integer(baseline[ids]),
  sex_residualized_cluster = as.integer(
    sex_assignments$sex_residualized_cluster
  ),
  clean_strip_residualized_cluster = as.integer(
    strip_assignments$clean_strip_residualized_cluster
  )
) |>
  dplyr::mutate(
    changed_after_sex_residualization =
      baseline_cluster != sex_residualized_cluster,
    changed_after_clean_strip_residualization =
      baseline_cluster != clean_strip_residualized_cluster,
    changed_in_either_sensitivity_analysis =
      changed_after_sex_residualization |
        changed_after_clean_strip_residualization,
    changed_in_both_sensitivity_analyses =
      changed_after_sex_residualization &
        changed_after_clean_strip_residualization,
    same_cluster_after_both_residualizations =
      sex_residualized_cluster == clean_strip_residualized_cluster,
    transition_pattern = paste0(
      baseline_cluster,
      " -> ",
      sex_residualized_cluster,
      " -> ",
      clean_strip_residualized_cluster
    )
  )

assignment_transitions_long <- assignment_transitions |>
  dplyr::select(
    participant_id,
    sex,
    baseline_cluster,
    sex_residualized_cluster,
    clean_strip_residualized_cluster
  ) |>
  tidyr::pivot_longer(
    cols = c(
      baseline_cluster,
      sex_residualized_cluster,
      clean_strip_residualized_cluster
    ),
    names_to = "analysis",
    values_to = "cluster"
  ) |>
  dplyr::mutate(
    analysis = factor(
      analysis,
      levels = c(
        "baseline_cluster",
        "sex_residualized_cluster",
        "clean_strip_residualized_cluster"
      ),
      labels = c(
        "Baseline",
        "Sex-residualized",
        "Clean-strip-residualized"
      )
    ),
    participant_order = match(participant_id, ids),
    analysis_order = as.integer(analysis),
    cluster = as.integer(cluster)
  ) |>
  dplyr::arrange(participant_order, analysis_order)


# ------------------------------------------------------------
# 4. Pairwise partition agreement
# ------------------------------------------------------------

partitions <- list(
  Baseline = assignment_transitions$baseline_cluster,
  Sex_residualized = assignment_transitions$sex_residualized_cluster,
  Clean_strip_residualized =
    assignment_transitions$clean_strip_residualized_cluster
)

partition_names <- names(partitions)

ari_matrix <- matrix(
  NA_real_,
  nrow = length(partitions),
  ncol = length(partitions),
  dimnames = list(partition_names, partition_names)
)

for (i in seq_along(partitions)) {
  for (j in seq_along(partitions)) {
    ari_matrix[i, j] <- mclust::adjustedRandIndex(
      partitions[[i]],
      partitions[[j]]
    )
  }
}

ari_matrix_wide <- as.data.frame(ari_matrix) |>
  tibble::rownames_to_column(var = "partition")

ari_matrix_long <- as.data.frame(as.table(ari_matrix)) |>
  tibble::as_tibble() |>
  dplyr::rename(
    partition_1 = Var1,
    partition_2 = Var2,
    ARI = Freq
  ) |>
  dplyr::mutate(
    partition_1 = as.character(partition_1),
    partition_2 = as.character(partition_2)
  )

pairwise_partition_summary <- dplyr::bind_rows(
  tibble::tibble(
    comparison = "Baseline vs Sex-residualized",
    ARI = mclust::adjustedRandIndex(
      partitions$Baseline,
      partitions$Sex_residualized
    ),
    same_aligned_label_n = sum(
      partitions$Baseline == partitions$Sex_residualized
    ),
    different_aligned_label_n = sum(
      partitions$Baseline != partitions$Sex_residualized
    )
  ),
  tibble::tibble(
    comparison = "Baseline vs Clean-strip-residualized",
    ARI = mclust::adjustedRandIndex(
      partitions$Baseline,
      partitions$Clean_strip_residualized
    ),
    same_aligned_label_n = sum(
      partitions$Baseline == partitions$Clean_strip_residualized
    ),
    different_aligned_label_n = sum(
      partitions$Baseline != partitions$Clean_strip_residualized
    )
  ),
  tibble::tibble(
    comparison = "Sex-residualized vs Clean-strip-residualized",
    ARI = mclust::adjustedRandIndex(
      partitions$Sex_residualized,
      partitions$Clean_strip_residualized
    ),
    same_aligned_label_n = sum(
      partitions$Sex_residualized ==
        partitions$Clean_strip_residualized
    ),
    different_aligned_label_n = sum(
      partitions$Sex_residualized !=
        partitions$Clean_strip_residualized
    )
  )
)

cluster_size_summary <- dplyr::bind_rows(
  lapply(
    names(partitions),
    function(partition_name) {
      tibble::tibble(
        partition = partition_name,
        cluster = 1:3,
        n = as.integer(
          table(factor(partitions[[partition_name]], levels = 1:3))
        )
      )
    }
  )
)


# ------------------------------------------------------------
# 5. Tidy contingency tables for later figures/tables
# ------------------------------------------------------------

make_contingency <- function(x, y, x_name, y_name) {
  as.data.frame(table(x, y)) |>
    tibble::as_tibble() |>
    dplyr::rename(
      cluster_1 = x,
      cluster_2 = y,
      n = Freq
    ) |>
    dplyr::mutate(
      comparison = paste(x_name, "vs", y_name),
      partition_1 = x_name,
      partition_2 = y_name,
      cluster_1 = as.integer(as.character(cluster_1)),
      cluster_2 = as.integer(as.character(cluster_2))
    ) |>
    dplyr::select(
      comparison,
      partition_1,
      partition_2,
      cluster_1,
      cluster_2,
      n
    )
}

partition_contingencies <- dplyr::bind_rows(
  make_contingency(
    factor(partitions$Baseline, levels = 1:3),
    factor(partitions$Sex_residualized, levels = 1:3),
    "Baseline",
    "Sex-residualized"
  ),
  make_contingency(
    factor(partitions$Baseline, levels = 1:3),
    factor(partitions$Clean_strip_residualized, levels = 1:3),
    "Baseline",
    "Clean-strip-residualized"
  ),
  make_contingency(
    factor(partitions$Sex_residualized, levels = 1:3),
    factor(partitions$Clean_strip_residualized, levels = 1:3),
    "Sex-residualized",
    "Clean-strip-residualized"
  )
)

transition_pattern_summary <- assignment_transitions |>
  dplyr::count(transition_pattern, name = "n") |>
  dplyr::arrange(dplyr::desc(n), transition_pattern)

participant_change_summary <- tibble::tibble(
  metric = c(
    "Participants changing after sex residualization",
    "Participants changing after clean-strip residualization",
    "Participants changing in either sensitivity analysis",
    "Participants changing in both sensitivity analyses",
    "Participants assigned to the same cluster after both residualizations"
  ),
  value = c(
    sum(assignment_transitions$changed_after_sex_residualization),
    sum(assignment_transitions$changed_after_clean_strip_residualization),
    sum(assignment_transitions$changed_in_either_sensitivity_analysis),
    sum(assignment_transitions$changed_in_both_sensitivity_analyses),
    sum(assignment_transitions$same_cluster_after_both_residualizations)
  )
)


# ------------------------------------------------------------
# 6. Cross-check against the source summary files
# ------------------------------------------------------------

sex_summary <- readr::read_csv(sex_summary_file, show_col_types = FALSE)
strip_summary <- readr::read_csv(strip_summary_file, show_col_types = FALSE)

sex_expected_ari <- as.numeric(sex_summary$ARI_vs_baseline[[1]])
strip_expected_ari <- as.numeric(
  strip_summary$value[
    strip_summary$metric == "ARI vs baseline"
  ][[1]]
)

if (
  abs(
    pairwise_partition_summary$ARI[
      pairwise_partition_summary$comparison ==
        "Baseline vs Sex-residualized"
    ] - sex_expected_ari
  ) > 1e-12
) {
  stop("The combined sex-residualized ARI does not match script 04.")
}

if (
  abs(
    pairwise_partition_summary$ARI[
      pairwise_partition_summary$comparison ==
        "Baseline vs Clean-strip-residualized"
    ] - strip_expected_ari
  ) > 1e-12
) {
  stop("The combined clean-strip-residualized ARI does not match script 05.")
}


# ------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  assignment_transitions,
  file.path(tables_dir, "assignment_transitions_participant_level.csv")
)

readr::write_csv(
  assignment_transitions_long,
  file.path(tables_dir, "assignment_transitions_long.csv")
)

readr::write_csv(
  ari_matrix_wide,
  file.path(tables_dir, "partition_ARI_matrix.csv")
)

readr::write_csv(
  ari_matrix_long,
  file.path(tables_dir, "partition_ARI_matrix_long.csv")
)

readr::write_csv(
  pairwise_partition_summary,
  file.path(tables_dir, "pairwise_partition_summary.csv")
)

readr::write_csv(
  cluster_size_summary,
  file.path(tables_dir, "partition_cluster_sizes.csv")
)

readr::write_csv(
  partition_contingencies,
  file.path(tables_dir, "partition_contingencies_long.csv")
)

readr::write_csv(
  transition_pattern_summary,
  file.path(tables_dir, "transition_pattern_summary.csv")
)

readr::write_csv(
  participant_change_summary,
  file.path(tables_dir, "participant_change_summary.csv")
)

assignment_transition_objects <- list(
  settings = list(
    seed = SEED,
    partitions = partition_names,
    cluster_labels_are_aligned_for_display = TRUE
  ),
  assignment_transitions = assignment_transitions,
  assignment_transitions_long = assignment_transitions_long,
  ari_matrix = ari_matrix,
  ari_matrix_wide = ari_matrix_wide,
  ari_matrix_long = ari_matrix_long,
  pairwise_partition_summary = pairwise_partition_summary,
  cluster_size_summary = cluster_size_summary,
  partition_contingencies = partition_contingencies,
  transition_pattern_summary = transition_pattern_summary,
  participant_change_summary = participant_change_summary
)

saveRDS(
  assignment_transition_objects,
  file.path(processed_dir, "assignment_transition_objects.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(tables_dir, "sessionInfo_assignment_transitions.txt")
)


# ------------------------------------------------------------
# 8. Console report
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ASSIGNMENT TRANSITIONS COMPLETE\n")
cat("========================================\n")
print(pairwise_partition_summary, n = Inf, width = Inf)
cat("\n")
print(participant_change_summary, n = Inf, width = Inf)
cat("========================================\n")
