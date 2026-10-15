# ============================================================
# validate_assignment_transition_outputs.R
# Numerical guardrails for R/06_assignment_transitions.R
# ============================================================

required_files <- c(
  file.path("tables", "assignment_transitions_participant_level.csv"),
  file.path("tables", "assignment_transitions_long.csv"),
  file.path("tables", "partition_ARI_matrix.csv"),
  file.path("tables", "pairwise_partition_summary.csv"),
  file.path("tables", "participant_change_summary.csv"),
  file.path(
    "data",
    "pediatric",
    "processed",
    "assignment_transition_objects.rds"
  )
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing output file(s). Run R/06_assignment_transitions.R first:\n",
    paste(missing_files, collapse = "\n")
  )
}

transitions <- read.csv(
  file.path("tables", "assignment_transitions_participant_level.csv"),
  check.names = FALSE
)

long_tbl <- read.csv(
  file.path("tables", "assignment_transitions_long.csv"),
  check.names = FALSE
)

pairwise <- read.csv(
  file.path("tables", "pairwise_partition_summary.csv"),
  check.names = FALSE
)

change_summary <- read.csv(
  file.path("tables", "participant_change_summary.csv"),
  check.names = FALSE
)

stopifnot(
  nrow(transitions) == 25L,
  length(unique(transitions$participant_id)) == 25L,
  nrow(long_tbl) == 75L,
  nrow(pairwise) == 3L,
  !anyNA(
    transitions[
      c(
        "baseline_cluster",
        "sex_residualized_cluster",
        "clean_strip_residualized_cluster"
      )
    ]
  )
)

sex_ari <- pairwise$ARI[
  pairwise$comparison == "Baseline vs Sex-residualized"
]
strip_ari <- pairwise$ARI[
  pairwise$comparison == "Baseline vs Clean-strip-residualized"
]

stopifnot(
  length(sex_ari) == 1L,
  length(strip_ari) == 1L,
  abs(sex_ari - 0.134) < 0.02,
  abs(strip_ari - 0.407) < 0.02
)

stopifnot(
  change_summary$value[
    change_summary$metric ==
      "Participants changing after sex residualization"
  ] == sum(transitions$changed_after_sex_residualization),
  change_summary$value[
    change_summary$metric ==
      "Participants changing after clean-strip residualization"
  ] == sum(transitions$changed_after_clean_strip_residualization)
)

obj <- readRDS(
  file.path(
    "data",
    "pediatric",
    "processed",
    "assignment_transition_objects.rds"
  )
)

ari <- obj$ari_matrix
stopifnot(
  identical(dim(ari), c(3L, 3L)),
  max(abs(ari - t(ari))) < 1e-12,
  max(abs(diag(ari) - 1)) < 1e-12
)

cat("Assignment-transition validation passed.\n")
