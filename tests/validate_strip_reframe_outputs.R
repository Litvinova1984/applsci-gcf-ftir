# ============================================================
# validate_strip_reframe_outputs.R
# Numerical guardrails for the cleaned strip-sensitivity module
# ============================================================
# Run after:
#   source("R/05_strip_reproducibility.R")
# from the repository root.
#
# Expected values come from the frozen pre-reframe public workflow.
# The cleaned module must preserve the independent clean-strip results
# while removing data-dependent interval/domain analyses and clinical
# reference-proximity branches.
# ============================================================

required_files <- c(
  file.path("tables", "strip_technical_reproducibility_summary.csv"),
  file.path("tables", "strip_reference_residualization_sensitivity.csv"),
  file.path("tables", "strip_reference_residualization_pairwise_ARI.csv"),
  file.path("tables", "strip_primary_sensitivity_summary.csv"),
  file.path("tables", "strip_primary_assignment_transition.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Run R/05_strip_reproducibility.R first. Missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

technical <- read.csv(
  file.path("tables", "strip_technical_reproducibility_summary.csv"),
  stringsAsFactors = FALSE
)
resid <- read.csv(
  file.path("tables", "strip_reference_residualization_sensitivity.csv"),
  stringsAsFactors = FALSE
)
pairwise <- read.csv(
  file.path("tables", "strip_reference_residualization_pairwise_ARI.csv"),
  stringsAsFactors = FALSE
)
primary <- read.csv(
  file.path("tables", "strip_primary_sensitivity_summary.csv"),
  stringsAsFactors = FALSE
)
transition <- read.csv(
  file.path("tables", "strip_primary_assignment_transition.csv"),
  stringsAsFactors = FALSE
)

metric_value <- function(tbl, name) {
  z <- tbl$value[tbl$metric == name]
  if (length(z) != 1L) {
    stop("Expected exactly one value for metric: ", name)
  }
  as.numeric(z)
}

expect_close <- function(observed, expected, tolerance, label) {
  if (!isTRUE(all.equal(observed, expected, tolerance = tolerance))) {
    stop(
      label,
      " mismatch. Observed = ",
      paste(observed, collapse = ", "),
      "; expected = ",
      paste(expected, collapse = ", ")
    )
  }
}

expect_close(
  metric_value(technical, "Minimum pairwise Pearson"),
  0.9999889915461593,
  1e-12,
  "Minimum pairwise Pearson"
)
expect_close(
  metric_value(technical, "Median pairwise Pearson"),
  0.9999955220471575,
  1e-12,
  "Median pairwise Pearson"
)
expect_close(
  metric_value(technical, "Minimum pairwise Spearman"),
  0.9999822469862572,
  1e-12,
  "Minimum pairwise Spearman"
)
expect_close(
  metric_value(technical, "Median pairwise Spearman"),
  0.9999900616807382,
  1e-12,
  "Median pairwise Spearman"
)
expect_close(
  metric_value(technical, "Integrated absorbance CV (%)"),
  0.8187663943672774,
  1e-10,
  "Integrated absorbance CV"
)
expect_close(
  metric_value(technical, "Median pointwise CV (%)"),
  0.8542025069922996,
  1e-10,
  "Median pointwise CV"
)
expect_close(
  metric_value(technical, "95th percentile pointwise CV (%)"),
  1.174046212992733,
  1e-10,
  "95th percentile pointwise CV"
)
expect_close(
  metric_value(technical, "Archived strip vs four-strip mean Pearson"),
  0.9999991825802927,
  1e-12,
  "Archived-vs-mean Pearson"
)
expect_close(
  metric_value(technical, "Archived strip vs four-strip mean Spearman"),
  0.9999957948683891,
  1e-12,
  "Archived-vs-mean Spearman"
)

if (nrow(resid) != 6L) {
  stop("Expected six clean-strip reference directions; found ", nrow(resid))
}
expect_close(
  resid$ARI_vs_baseline,
  rep(0.4071146245059289, 6L),
  1e-12,
  "Reference-specific ARI"
)

if (!all(apply(resid[, c("cluster_1_n", "cluster_2_n", "cluster_3_n")], 1, function(x) {
  identical(sort(as.integer(x), decreasing = TRUE), c(11L, 8L, 6L))
}))) {
  stop("Residualized cluster sizes are not 11 / 8 / 6 for every reference.")
}

expect_close(
  pairwise$ARI,
  rep(1, nrow(pairwise)),
  1e-12,
  "Pairwise ARI among clean-strip residualized partitions"
)

expect_close(
  metric_value(primary, "ARI vs baseline"),
  0.4071146245059289,
  1e-12,
  "Primary ARI"
)
expect_close(
  metric_value(primary, "Three-PC cumulative variance (%)"),
  89.8237650611702,
  1e-10,
  "Primary three-PC cumulative variance"
)
expect_close(
  metric_value(primary, "Average silhouette"),
  0.37126592298544175,
  1e-10,
  "Primary average silhouette"
)
expect_close(
  metric_value(primary, "Negative silhouette count"),
  2,
  0,
  "Primary negative silhouette count"
)

if (nrow(transition) != 25L) {
  stop("Expected 25 participant-level transition rows; found ", nrow(transition))
}
if (!all(transition$baseline_cluster %in% 1:3)) {
  stop("Invalid baseline cluster labels in transition table.")
}
if (!all(transition$clean_strip_residualized_cluster %in% 1:3)) {
  stop("Invalid residualized cluster labels in transition table.")
}

cat("Strip reframe numerical guardrails: PASS\n")
