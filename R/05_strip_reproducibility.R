# ============================================================
# 05_strip_reproducibility.R
# Pediatric GCF FTIR exploratory reframe
# Clean-strip reproducibility and substrate-direction sensitivity
# ============================================================
#
# This script retains only technical references and spectral directions
# defined independently of the data-derived pediatric cluster labels:
#
#   1. reproducibility of four independent unused strips;
#   2. agreement of those strips with the archived clean-strip spectrum;
#   3. projection of clean-strip references as supplementary observations
#      into the fixed pediatric PCA space;
#   4. sensitivity of the fixed pediatric PCA/HCA procedure to removal of
#      a clean-strip-aligned spectral direction.
#
# The four DPT files are already strip-level averaged spectra (three
# spatial measurement locations per strip). Point-level measurements are
# unavailable and are not treated as independent observations.
#
# Clean strips are technical references, not biological controls. Removing
# a clean-strip-aligned direction is a sensitivity analysis, not a validated
# correction that separates technical and biological signal.
#
# No data-dependent spectral intervals, cluster-wise inferential testing,
# supplementary pediatric reference proximity, or interval chosen from the
# cluster solution is used here.
# ============================================================

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}
source(file.path("R", "00_utils.R"))

required_pkgs <- c(
  "FactoMineR",
  "mclust",
  "cluster",
  "tibble",
  "dplyr",
  "readr",
  "ggplot2"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

processed_dir <- file.path("data", "pediatric", "processed")
technical_dir <- file.path("data", "technical_controls")
tables_dir <- "tables"
figures_dir <- "figures"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(
  file.path(processed_dir, "pediatric_X_snv.rds"),
  file.path(processed_dir, "pediatric_wavenumbers.rds"),
  file.path(processed_dir, "pediatric_source_spectra_long.rds"),
  file.path(processed_dir, "pediatric_clusters_baseline.rds")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Required processed input(s) are missing:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

X_snv <- as.matrix(
  readRDS(file.path(processed_dir, "pediatric_X_snv.rds"))
)
wn <- as.numeric(
  readRDS(file.path(processed_dir, "pediatric_wavenumbers.rds"))
)
spectra_long <- readRDS(
  file.path(processed_dir, "pediatric_source_spectra_long.rds")
)
cluster_raw <- readRDS(
  file.path(processed_dir, "pediatric_clusters_baseline.rds")
)

if (!identical(dim(X_snv), c(25L, 793L))) {
  stop("Expected X_snv = 25 x 793.")
}
if (length(wn) != ncol(X_snv)) {
  stop("Wavenumber vector does not match X_snv.")
}

ids <- rownames(X_snv)
if (is.null(ids)) {
  stop("X_snv must have participant IDs as row names.")
}

if (is.null(names(cluster_raw))) {
  cluster_baseline <- as.integer(cluster_raw)
  names(cluster_baseline) <- ids
} else {
  cluster_raw <- cluster_raw[ids]
  cluster_baseline <- if (is.factor(cluster_raw)) {
    as.integer(as.character(cluster_raw))
  } else {
    as.integer(cluster_raw)
  }
  names(cluster_baseline) <- ids
}

if (!identical(as.integer(table(cluster_baseline)), c(15L, 4L, 6L))) {
  stop("Expected baseline cluster sizes 15 / 4 / 6.")
}

N_PC_BASELINE <- 3L
K_BASELINE <- 3L
GRID_TOLERANCE_CM1 <- 1e-4

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

ward_k3 <- function(scores) {
  stats::cutree(
    stats::hclust(
      stats::dist(scores, method = "euclidean"),
      method = "ward.D2"
    ),
    k = K_BASELINE
  )
}

align_cluster_labels <- function(reference, target) {
  reference <- as.integer(reference)
  target <- as.integer(target)

  if (!identical(sort(unique(reference)), seq_len(K_BASELINE))) {
    stop("Reference labels must be 1, 2, and 3.")
  }
  if (!identical(sort(unique(target)), seq_len(K_BASELINE))) {
    stop("Target labels must be 1, 2, and 3.")
  }

  permutations <- rbind(
    c(1L, 2L, 3L),
    c(1L, 3L, 2L),
    c(2L, 1L, 3L),
    c(2L, 3L, 1L),
    c(3L, 1L, 2L),
    c(3L, 2L, 1L)
  )

  agreements <- apply(
    permutations,
    1,
    function(p) sum(p[target] == reference)
  )

  best <- permutations[which.max(agreements), ]
  aligned <- as.integer(best[target])
  names(aligned) <- names(target)
  aligned
}

trapz_base <- function(x, y) {
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

integrated_two_windows <- function(x, y) {
  low <- x >= analysis_windows$low[1] & x <= analysis_windows$low[2]
  high <- x >= analysis_windows$high[1] & x <= analysis_windows$high[2]
  trapz_base(x[low], y[low]) + trapz_base(x[high], y[high])
}

cv_percent <- function(x) {
  100 * stats::sd(x) / abs(mean(x))
}

read_dpt_spectrum <- function(path) {
  z <- utils::read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("wavenumber", "absorbance")
  )

  if (ncol(z) != 2L) {
    stop("Expected exactly two comma-separated columns in ", path)
  }

  out <- tibble::tibble(
    wavenumber = suppressWarnings(as.numeric(z$wavenumber)),
    absorbance = suppressWarnings(as.numeric(z$absorbance))
  ) |>
    dplyr::filter(is.finite(wavenumber), is.finite(absorbance)) |>
    dplyr::arrange(dplyr::desc(wavenumber))

  if (nrow(out) < 100L) {
    stop("Too few spectral points parsed from ", path, ": ", nrow(out))
  }

  out
}

extract_archived_clean_strip <- function() {
  z <- spectra_long |>
    dplyr::filter(
      sample_id == "etalon",
      retain_analysis_windows(wavenumber)
    ) |>
    dplyr::arrange(dplyr::desc(wavenumber))

  if (nrow(z) != length(wn)) {
    stop("Unexpected retained grid for archived clean strip.")
  }
  if (!isTRUE(all.equal(as.numeric(z$wavenumber), wn, tolerance = 1e-12))) {
    stop("Wavenumber-grid mismatch for archived clean strip.")
  }

  as.numeric(z$absorbance)
}

remove_reference_component <- function(X, reference_vector) {
  reference_vector <- as.numeric(reference_vector)

  if (length(reference_vector) != ncol(X)) {
    stop("Reference vector length does not match spectral matrix.")
  }

  e <- reference_vector - mean(reference_vector)
  e_norm <- sqrt(sum(e * e))

  if (!is.finite(e_norm) || e_norm <= 0) {
    stop("Reference vector has zero or non-finite centered norm.")
  }

  e_unit <- e / e_norm
  X_centered <- sweep(X, 1, rowMeans(X), "-")
  coefficients <- as.numeric(X_centered %*% e_unit)
  X_residual <- X_centered - tcrossprod(coefficients, e_unit)

  rownames(X_residual) <- rownames(X)
  colnames(X_residual) <- colnames(X)
  names(coefficients) <- rownames(X)

  list(
    matrix = X_residual,
    coefficients = coefficients,
    unit_reference = e_unit
  )
}

# ============================================================
# 1. Archived and four independent clean-strip references
# ============================================================

archived_strip_raw <- extract_archived_clean_strip()
archived_strip_snv <- snv(archived_strip_raw)

strip_files <- file.path(
  technical_dir,
  sprintf("etalon_%02d.dpt", 1:4)
)
if (!all(file.exists(strip_files))) {
  stop("Expected: ", paste(strip_files, collapse = ", "))
}

strip_list <- lapply(strip_files, read_dpt_spectrum)
names(strip_list) <- sprintf("strip_%02d", 1:4)
strip_list <- lapply(
  strip_list,
  function(z) {
    z |>
      dplyr::filter(retain_analysis_windows(wavenumber)) |>
      dplyr::arrange(dplyr::desc(wavenumber))
  }
)

lens <- vapply(strip_list, nrow, integer(1))
if (length(unique(lens)) != 1L) {
  stop("Four strips have different retained grid lengths.")
}

strip_wn <- strip_list[[1]]$wavenumber
for (i in seq_along(strip_list)) {
  if (!isTRUE(all.equal(strip_list[[i]]$wavenumber, strip_wn, tolerance = 1e-12))) {
    stop("Four strips have different wavenumber grids.")
  }
}

if (length(strip_wn) != length(wn)) {
  stop("DPT and archived pediatric grids have different lengths.")
}

max_grid_difference <- max(abs(strip_wn - wn))
if (!is.finite(max_grid_difference) || max_grid_difference > GRID_TOLERANCE_CM1) {
  stop(
    "DPT grid differs from the archived pediatric grid by more than ",
    GRID_TOLERANCE_CM1,
    " cm^-1. Maximum difference = ",
    signif(max_grid_difference, 6),
    " cm^-1."
  )
}

X_strips <- do.call(
  rbind,
  lapply(strip_list, function(z) z$absorbance)
)
rownames(X_strips) <- names(strip_list)
colnames(X_strips) <- colnames(X_snv)

X_strips_snv <- snv_matrix(X_strips)
consensus_four_strip_snv <- colMeans(X_strips_snv)

# ============================================================
# 2. Technical reproducibility of unused strips
# ============================================================

pearson_matrix <- stats::cor(t(X_strips), method = "pearson")
spearman_matrix <- stats::cor(t(X_strips), method = "spearman")
pairs <- which(upper.tri(pearson_matrix), arr.ind = TRUE)

pairwise_correlations <- tibble::tibble(
  strip_1 = rownames(pearson_matrix)[pairs[, 1]],
  strip_2 = colnames(pearson_matrix)[pairs[, 2]],
  pearson_r = pearson_matrix[pairs],
  spearman_rho = spearman_matrix[pairs]
)

integrated_absorbance <- tibble::tibble(
  strip = rownames(X_strips),
  integrated_absorbance = apply(
    X_strips,
    1,
    function(y) integrated_two_windows(strip_wn, y)
  )
)
integrated_cv <- cv_percent(integrated_absorbance$integrated_absorbance)

pointwise_mean <- colMeans(X_strips)
pointwise_sd <- apply(X_strips, 2, stats::sd)
pointwise_cv <- 100 * pointwise_sd / abs(pointwise_mean)
pointwise_cv[!is.finite(pointwise_cv)] <- NA_real_

pointwise_variability <- tibble::tibble(
  wavenumber = strip_wn,
  mean_absorbance = pointwise_mean,
  sd_absorbance = pointwise_sd,
  cv_percent = pointwise_cv
)

pointwise_summary <- tibble::tibble(
  metric = c(
    "Median pointwise CV (%)",
    "Mean pointwise CV (%)",
    "95th percentile pointwise CV (%)",
    "Maximum pointwise CV (%)"
  ),
  value = c(
    stats::median(pointwise_cv, na.rm = TRUE),
    mean(pointwise_cv, na.rm = TRUE),
    stats::quantile(pointwise_cv, 0.95, na.rm = TRUE, names = FALSE),
    max(pointwise_cv, na.rm = TRUE)
  )
)

archived_vs_four <- tibble::tibble(
  strip = rownames(X_strips),
  pearson_r = apply(
    X_strips,
    1,
    stats::cor,
    y = archived_strip_raw,
    method = "pearson"
  ),
  spearman_rho = apply(
    X_strips,
    1,
    stats::cor,
    y = archived_strip_raw,
    method = "spearman"
  )
)

four_strip_raw_mean <- colMeans(X_strips)
archived_vs_mean <- tibble::tibble(
  pearson_r = stats::cor(
    archived_strip_raw,
    four_strip_raw_mean,
    method = "pearson"
  ),
  spearman_rho = stats::cor(
    archived_strip_raw,
    four_strip_raw_mean,
    method = "spearman"
  )
)

technical_summary <- tibble::tibble(
  metric = c(
    "Minimum pairwise Pearson",
    "Median pairwise Pearson",
    "Minimum pairwise Spearman",
    "Median pairwise Spearman",
    "Integrated absorbance CV (%)",
    "Median pointwise CV (%)",
    "95th percentile pointwise CV (%)",
    "Archived strip vs four-strip mean Pearson",
    "Archived strip vs four-strip mean Spearman",
    "Maximum DPT-vs-archived grid difference (cm^-1)"
  ),
  value = c(
    min(pairwise_correlations$pearson_r),
    stats::median(pairwise_correlations$pearson_r),
    min(pairwise_correlations$spearman_rho),
    stats::median(pairwise_correlations$spearman_rho),
    integrated_cv,
    pointwise_summary$value[
      pointwise_summary$metric == "Median pointwise CV (%)"
    ],
    pointwise_summary$value[
      pointwise_summary$metric == "95th percentile pointwise CV (%)"
    ],
    archived_vs_mean$pearson_r,
    archived_vs_mean$spearman_rho,
    max_grid_difference
  )
)

# ============================================================
# 3. Supplementary projection of clean-strip references
# ============================================================
#
# The pediatric observations alone define the PCA axes. Clean-strip
# references are supplementary observations and do not affect the model.

reference_vectors_snv <- rbind(
  archived_clean_strip = archived_strip_snv,
  X_strips_snv,
  mean_four_strips = consensus_four_strip_snv
)
colnames(reference_vectors_snv) <- colnames(X_snv)

X_projection <- rbind(X_snv, reference_vectors_snv)
ind_sup <- (nrow(X_snv) + 1L):nrow(X_projection)

pca_clean_projection <- FactoMineR::PCA(
  X_projection,
  scale.unit = FALSE,
  ncp = N_PC_BASELINE,
  ind.sup = ind_sup,
  graph = FALSE
)

active_projection <- as.matrix(
  pca_clean_projection$ind$coord[, 1:N_PC_BASELINE, drop = FALSE]
)
reference_projection <- as.matrix(
  pca_clean_projection$ind.sup$coord[, 1:N_PC_BASELINE, drop = FALSE]
)
rownames(active_projection) <- ids
rownames(reference_projection) <- rownames(reference_vectors_snv)

clean_reference_projection <- tibble::tibble(
  reference = rownames(reference_projection),
  PC1 = reference_projection[, 1],
  PC2 = reference_projection[, 2],
  PC3 = reference_projection[, 3]
)

# ============================================================
# 4. Sensitivity to removal of a clean-strip-aligned direction
# ============================================================
#
# The same fixed procedure is repeated for every independently measured
# clean-strip reference: three PCs, Euclidean distance, Ward.D2, k = 3.
# Label alignment is used only to display participant transitions; ARI is
# label invariant.

residualization_objects <- lapply(
  seq_len(nrow(reference_vectors_snv)),
  function(i) {
    reference_name <- rownames(reference_vectors_snv)[i]

    residualized <- remove_reference_component(
      X_snv,
      reference_vectors_snv[i, , drop = TRUE]
    )

    pca_residualized <- FactoMineR::PCA(
      residualized$matrix,
      scale.unit = FALSE,
      ncp = N_PC_BASELINE,
      graph = FALSE
    )

    scores_residualized <- as.matrix(
      pca_residualized$ind$coord[
        ,
        1:N_PC_BASELINE,
        drop = FALSE
      ]
    )
    rownames(scores_residualized) <- ids

    cluster_unaligned <- ward_k3(scores_residualized)
    names(cluster_unaligned) <- ids
    cluster_aligned <- align_cluster_labels(
      cluster_baseline,
      cluster_unaligned
    )
    names(cluster_aligned) <- ids

    silhouette_object <- cluster::silhouette(
      cluster_unaligned,
      stats::dist(scores_residualized, method = "euclidean")
    )

    list(
      reference = reference_name,
      residual_matrix = residualized$matrix,
      projection_coefficients = residualized$coefficients,
      pca = pca_residualized,
      scores = scores_residualized,
      cluster_unaligned = cluster_unaligned,
      cluster_aligned = cluster_aligned,
      ARI_vs_baseline = mclust::adjustedRandIndex(
        cluster_baseline,
        cluster_unaligned
      ),
      changed_participants = sum(cluster_aligned != cluster_baseline),
      cumulative_variance_3PC_percent = sum(
        pca_residualized$eig[1:N_PC_BASELINE, 2]
      ),
      average_silhouette = mean(silhouette_object[, "sil_width"]),
      negative_silhouette_count = sum(
        silhouette_object[, "sil_width"] < 0
      )
    )
  }
)
names(residualization_objects) <- rownames(reference_vectors_snv)

reference_residualization_summary <- dplyr::bind_rows(
  lapply(
    residualization_objects,
    function(obj) {
      sizes <- as.integer(
        table(factor(obj$cluster_aligned, levels = 1:K_BASELINE))
      )

      tibble::tibble(
        reference = obj$reference,
        ARI_vs_baseline = obj$ARI_vs_baseline,
        changed_participants = obj$changed_participants,
        cumulative_variance_3PC_percent =
          obj$cumulative_variance_3PC_percent,
        average_silhouette = obj$average_silhouette,
        negative_silhouette_count = obj$negative_silhouette_count,
        cluster_1_n = sizes[1],
        cluster_2_n = sizes[2],
        cluster_3_n = sizes[3]
      )
    }
  )
)

reference_names <- names(residualization_objects)
reference_residualization_pairwise_ARI <- dplyr::bind_rows(
  lapply(
    seq_along(reference_names),
    function(i) {
      if (i == length(reference_names)) {
        return(NULL)
      }

      dplyr::bind_rows(
        lapply(
          (i + 1L):length(reference_names),
          function(j) {
            tibble::tibble(
              reference_1 = reference_names[i],
              reference_2 = reference_names[j],
              ARI = mclust::adjustedRandIndex(
                residualization_objects[[reference_names[i]]]$cluster_unaligned,
                residualization_objects[[reference_names[j]]]$cluster_unaligned
              )
            )
          }
        )
      )
    }
  )
)

primary_object <- residualization_objects[["mean_four_strips"]]

primary_assignment_transition <- tibble::tibble(
  participant_id = ids,
  baseline_cluster = as.integer(cluster_baseline),
  clean_strip_residualized_cluster = as.integer(
    primary_object$cluster_aligned[ids]
  ),
  assignment_changed = baseline_cluster != clean_strip_residualized_cluster,
  clean_strip_projection_coefficient = as.numeric(
    primary_object$projection_coefficients[ids]
  ),
  spearman_to_mean_clean_strip = apply(
    X_snv,
    1,
    stats::cor,
    y = consensus_four_strip_snv,
    method = "spearman"
  )
)

primary_transition_contingency <- as.data.frame.matrix(
  table(
    Baseline = primary_assignment_transition$baseline_cluster,
    CleanStripResidualized =
      primary_assignment_transition$clean_strip_residualized_cluster
  )
)
primary_transition_contingency <- tibble::rownames_to_column(
  primary_transition_contingency,
  var = "baseline_cluster"
)

primary_sensitivity_summary <- tibble::tibble(
  metric = c(
    "ARI vs baseline",
    "Participants changing aligned assignment",
    "Three-PC cumulative variance (%)",
    "Average silhouette",
    "Negative silhouette count"
  ),
  value = c(
    primary_object$ARI_vs_baseline,
    primary_object$changed_participants,
    primary_object$cumulative_variance_3PC_percent,
    primary_object$average_silhouette,
    primary_object$negative_silhouette_count
  )
)

# ============================================================
# 5. Diagnostic figures for later manuscript assembly
# ============================================================

clean_strip_overlay <- dplyr::bind_rows(
  tibble::tibble(
    wavenumber = wn,
    absorbance = archived_strip_raw,
    reference = "archived_clean_strip"
  ),
  tibble::tibble(
    wavenumber = rep(strip_wn, times = nrow(X_strips)),
    absorbance = as.numeric(t(X_strips)),
    reference = rep(rownames(X_strips), each = length(strip_wn))
  )
)

p_overlay <- ggplot2::ggplot(
  clean_strip_overlay,
  ggplot2::aes(wavenumber, absorbance, group = reference)
) +
  ggplot2::geom_line(linewidth = 0.5, alpha = 0.8) +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    title = "Archived and independent unused-strip spectra",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "Absorbance"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "technical_controls_clean_strip_overlay.png"),
  p_overlay,
  width = 10,
  height = 5.5,
  dpi = 300
)

p_cv <- ggplot2::ggplot(
  pointwise_variability,
  ggplot2::aes(wavenumber, cv_percent)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::scale_x_reverse() +
  ggplot2::labs(
    title = "Pointwise variability across four independent unused strips",
    x = expression("Wavenumber (cm"^{-1} * ")"),
    y = "Pointwise CV (%)"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "technical_controls_pointwise_CV.png"),
  p_cv,
  width = 10,
  height = 5.5,
  dpi = 300
)

projection_plot_data <- dplyr::bind_rows(
  tibble::tibble(
    id = ids,
    type = "pediatric_measurement",
    PC1 = active_projection[, 1],
    PC2 = active_projection[, 2]
  ),
  tibble::tibble(
    id = rownames(reference_projection),
    type = "clean_strip_reference",
    PC1 = reference_projection[, 1],
    PC2 = reference_projection[, 2]
  )
)

p_projection <- ggplot2::ggplot(
  projection_plot_data,
  ggplot2::aes(PC1, PC2, shape = type)
) +
  ggplot2::geom_point(size = 2.4, alpha = 0.8) +
  ggplot2::labs(
    title = "Clean-strip references projected into the pediatric PCA space",
    x = "PC1 score",
    y = "PC2 score",
    shape = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "technical_controls_PCA_projection.png"),
  p_projection,
  width = 7.5,
  height = 6,
  dpi = 300
)

# ============================================================
# 6. Save outputs
# ============================================================

readr::write_csv(
  pairwise_correlations,
  file.path(tables_dir, "strip_four_controls_pairwise_correlations.csv")
)
readr::write_csv(
  integrated_absorbance,
  file.path(tables_dir, "strip_four_controls_integrated_absorbance.csv")
)
readr::write_csv(
  pointwise_variability,
  file.path(tables_dir, "strip_four_controls_pointwise_variability.csv")
)
readr::write_csv(
  pointwise_summary,
  file.path(tables_dir, "strip_four_controls_pointwise_CV_summary.csv")
)
readr::write_csv(
  archived_vs_four,
  file.path(tables_dir, "strip_archived_vs_four_controls.csv")
)
readr::write_csv(
  technical_summary,
  file.path(tables_dir, "strip_technical_reproducibility_summary.csv")
)
readr::write_csv(
  clean_reference_projection,
  file.path(tables_dir, "strip_clean_reference_PCA_projection.csv")
)
readr::write_csv(
  reference_residualization_summary,
  file.path(tables_dir, "strip_reference_residualization_sensitivity.csv")
)
readr::write_csv(
  reference_residualization_pairwise_ARI,
  file.path(tables_dir, "strip_reference_residualization_pairwise_ARI.csv")
)
readr::write_csv(
  primary_assignment_transition,
  file.path(tables_dir, "strip_primary_assignment_transition.csv")
)
readr::write_csv(
  primary_transition_contingency,
  file.path(tables_dir, "strip_primary_transition_contingency.csv")
)
readr::write_csv(
  primary_sensitivity_summary,
  file.path(tables_dir, "strip_primary_sensitivity_summary.csv")
)

saveRDS(
  list(
    settings = list(
      baseline_n_pc = N_PC_BASELINE,
      baseline_k = K_BASELINE,
      primary_reference = "mean_four_strips",
      grid_tolerance_cm1 = GRID_TOLERANCE_CM1
    ),
    archived_strip_raw = archived_strip_raw,
    archived_strip_snv = archived_strip_snv,
    X_strips = X_strips,
    X_strips_snv = X_strips_snv,
    consensus_four_strip_snv = consensus_four_strip_snv,
    pairwise_correlations = pairwise_correlations,
    integrated_absorbance = integrated_absorbance,
    pointwise_variability = pointwise_variability,
    pointwise_summary = pointwise_summary,
    archived_vs_four = archived_vs_four,
    technical_summary = technical_summary,
    pca_clean_projection = pca_clean_projection,
    clean_reference_projection = clean_reference_projection,
    residualization_objects = residualization_objects,
    reference_residualization_summary = reference_residualization_summary,
    reference_residualization_pairwise_ARI =
      reference_residualization_pairwise_ARI,
    primary_assignment_transition = primary_assignment_transition,
    primary_transition_contingency = primary_transition_contingency,
    primary_sensitivity_summary = primary_sensitivity_summary
  ),
  file.path(processed_dir, "strip_reproducibility_objects.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(tables_dir, "sessionInfo_strip_reproducibility.txt")
)

cat("\n========================================\n")
cat("CLEAN-STRIP SENSITIVITY ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Maximum grid difference: ", signif(max_grid_difference, 6), " cm^-1\n", sep = "")
cat(
  "Pairwise Pearson range: ",
  paste(signif(range(pairwise_correlations$pearson_r), 8), collapse = " to "),
  "\n",
  sep = ""
)
cat(
  "Pairwise Spearman range: ",
  paste(signif(range(pairwise_correlations$spearman_rho), 8), collapse = " to "),
  "\n",
  sep = ""
)
cat("Integrated absorbance CV: ", round(integrated_cv, 3), "%\n", sep = "")
cat(
  "Median pointwise CV: ",
  round(stats::median(pointwise_cv, na.rm = TRUE), 3),
  "%\n",
  sep = ""
)
cat(
  "95th percentile pointwise CV: ",
  round(stats::quantile(pointwise_cv, 0.95, na.rm = TRUE), 3),
  "%\n",
  sep = ""
)
cat("\nResidualization using each clean-strip reference:\n")
print(reference_residualization_summary, n = Inf, width = Inf)
cat("\nPrimary mean-four-strip sensitivity:\n")
print(primary_sensitivity_summary, n = Inf, width = Inf)
cat("\nAssignment-transition contingency table:\n")
print(tibble::as_tibble(primary_transition_contingency), n = Inf, width = Inf)
cat("========================================\n")
