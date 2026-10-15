# ============================================================
# validate_pca_profiles_outputs.R
# Structural and numerical guardrails for R/02_pca_clustering.R
# ============================================================
# Run after:
#   source("R/02_pca_clustering.R")
# from the repository root.
# ============================================================

required_files <- c(
  file.path("tables", "PCA_eigenvalues_full.csv"),
  file.path("tables", "PCA_component_retention.csv"),
  file.path("tables", "PCA_scores_PC1_PC3.csv"),
  file.path("tables", "PCA_variable_coordinates_PC1_PC3.csv"),
  file.path("tables", "PCA_eigenvector_loadings_PC1_PC3.csv"),
  file.path("tables", "baseline_cluster_membership.csv"),
  file.path("tables", "baseline_cluster_sizes.csv"),
  file.path(
    "data",
    "pediatric",
    "processed",
    "pediatric_PCA_profiles_PC1_PC3.rds"
  )
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing output file(s). Run R/02_pca_clustering.R first:\n",
    paste(missing_files, collapse = "\n")
  )
}

eigenvalues <- read.csv(
  file.path("tables", "PCA_eigenvalues_full.csv"),
  check.names = FALSE
)

scores <- read.csv(
  file.path("tables", "PCA_scores_PC1_PC3.csv"),
  check.names = FALSE
)

varcoord <- read.csv(
  file.path("tables", "PCA_variable_coordinates_PC1_PC3.csv"),
  check.names = FALSE
)

loadings <- read.csv(
  file.path("tables", "PCA_eigenvector_loadings_PC1_PC3.csv"),
  check.names = FALSE
)

membership <- read.csv(
  file.path("tables", "baseline_cluster_membership.csv"),
  check.names = FALSE
)

sizes <- read.csv(
  file.path("tables", "baseline_cluster_sizes.csv"),
  check.names = FALSE
)

expected_profile_columns <- c(
  "wavenumber",
  "variable",
  "PC1",
  "PC2",
  "PC3"
)

stopifnot(
  nrow(scores) == 25L,
  length(unique(scores$sample_id)) == 25L,
  identical(names(scores), c("sample_id", "PC1", "PC2", "PC3")),
  nrow(varcoord) == 793L,
  nrow(loadings) == 793L,
  identical(names(varcoord), expected_profile_columns),
  identical(names(loadings), expected_profile_columns),
  nrow(membership) == 25L,
  length(unique(membership$sample_id)) == 25L,
  sum(sizes$n) == 25L,
  identical(sort(as.integer(sizes$n)), c(4L, 6L, 15L)),
  all(diff(varcoord$wavenumber) <= 0),
  all(diff(loadings$wavenumber) <= 0),
  identical(varcoord$variable, loadings$variable),
  max(abs(varcoord$wavenumber - loadings$wavenumber)) < 1e-12,
  all(is.finite(as.matrix(varcoord[c("PC1", "PC2", "PC3")]))),
  all(is.finite(as.matrix(loadings[c("PC1", "PC2", "PC3")]))),
  abs(eigenvalues$cumulative_percent[3] - 89.281825) < 0.05
)

V <- as.matrix(
  loadings[c("PC1", "PC2", "PC3")]
)

# Unit-norm eigenvectors should be orthonormal up to numerical error.
V_cross <- crossprod(V)
stopifnot(
  max(abs(diag(V_cross) - 1)) < 1e-8,
  max(abs(V_cross - diag(3))) < 1e-8
)

profiles <- readRDS(
  file.path(
    "data",
    "pediatric",
    "processed",
    "pediatric_PCA_profiles_PC1_PC3.rds"
  )
)

stopifnot(
  identical(dim(profiles$variable_coordinates), c(793L, 3L)),
  identical(dim(profiles$eigenvector_loadings), c(793L, 3L)),
  grepl("Raw mathematical PCA signs", profiles$sign_convention, fixed = TRUE)
)

cat("PCA profile validation passed.\n")
