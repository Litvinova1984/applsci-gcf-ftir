# ============================================================
# 05_strip_reproducibility.R
# Pediatric GCF FTIR
# Archived clean-strip sensitivity + four-strip reproducibility
# ============================================================
#
# PART A reproduces the technical-reference analyses based on the
# single archived clean unused strip included with the pediatric
# dataset. PART B evaluates four different physical unused strips
# from the same batch/instrumental series. PART C tests whether the
# substrate-related conclusions depend on the choice of a single
# clean-strip reference and whether multivariate structure remains
# after removal of a consensus four-strip component.
#
# The four DPT files are already averaged spectra (three spatial
# measurement points per strip). Point-level spectra are unavailable
# and are not treated as independent observations.
#
# The archived strip and the four additional strips are technical
# references, not biological controls.
#
# Reference distances are reproduced in the first five PCA dimensions
# because that was the specification used for that descriptive analysis.
# Baseline pediatric clustering itself remains the fixed 3-PC,
# Euclidean, Ward.D2, k = 3 solution established in R/02.
#
# Strip-removal and D1-region-exclusion analyses recompute PCA and
# fixed-k clustering and compare the new partition with baseline by ARI.
# Additional k = 2-6 indices after consensus-strip removal characterize
# residual structure; they are not treated as independent validation.
# ============================================================

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}
source(file.path("R", "00_utils.R"))

required_pkgs <- c("FactoMineR", "mclust", "cluster", "tibble", "dplyr", "readr", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))

processed_dir <- file.path("data", "pediatric", "processed")
technical_dir <- file.path("data", "technical_controls")
tables_dir <- "tables"
figures_dir <- "figures"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

X_snv <- as.matrix(readRDS(file.path(processed_dir, "pediatric_X_snv.rds")))
wn <- as.numeric(readRDS(file.path(processed_dir, "pediatric_wavenumbers.rds")))
spectra_long <- readRDS(file.path(processed_dir, "pediatric_source_spectra_long.rds"))
cluster_raw <- readRDS(file.path(processed_dir, "pediatric_clusters_baseline.rds"))

if (!identical(dim(X_snv), c(25L, 793L))) stop("Expected X_snv = 25 x 793.")
if (length(wn) != ncol(X_snv)) stop("Wavenumber vector does not match X_snv.")
ids <- rownames(X_snv)
if (is.null(ids)) stop("X_snv must have participant IDs as row names.")

if (is.null(names(cluster_raw))) {
  cluster_baseline <- as.integer(cluster_raw)
  names(cluster_baseline) <- ids
} else {
  cluster_raw <- cluster_raw[ids]
  cluster_baseline <- if (is.factor(cluster_raw)) as.integer(as.character(cluster_raw)) else as.integer(cluster_raw)
  names(cluster_baseline) <- ids
}
if (!identical(as.integer(table(cluster_baseline)), c(15L, 4L, 6L))) stop("Expected baseline sizes 15 / 4 / 6.")

N_PC_BASELINE <- 3L
N_PC_REFERENCE <- 5L
K_BASELINE <- 3L

ward_k3 <- function(scores) {
  stats::cutree(
    stats::hclust(stats::dist(scores, method = "euclidean"), method = "ward.D2"),
    k = K_BASELINE
  )
}

trapz_base <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

integrated_two_windows <- function(x, y) {
  low <- x >= 870 & x <= 1800
  high <- x >= 2800 & x <= 3400
  trapz_base(x[low], y[low]) + trapz_base(x[high], y[high])
}


cv_percent <- function(x) 100 * stats::sd(x) / abs(mean(x))


remove_reference_component <- function(X, reference_vector) {

  reference_vector <- as.numeric(reference_vector)

  if (length(reference_vector) != ncol(X)) {
    stop("Reference vector length does not match spectral matrix.")
  }

  e_ref <- reference_vector -
    mean(reference_vector)

  denom <- sum(e_ref * e_ref)

  if (!is.finite(denom) || denom <= 0) {
    stop("Reference vector has zero or non-finite centered variance.")
  }

  X_ctr <- t(
    apply(
      X,
      1,
      function(x) x - mean(x)
    )
  )

  beta_ref <- apply(
    X_ctr,
    1,
    function(x) {
      sum(x * e_ref) / denom
    }
  )

  X_res <- X_ctr -
    beta_ref %*%
    t(e_ref)

  rownames(X_res) <- rownames(X)
  colnames(X_res) <- colnames(X)

  X_res
}


calinski_harabasz <- function(scores, cl) {

  scores <- as.matrix(scores)
  cl <- factor(cl)

  n <- nrow(scores)
  k <- nlevels(cl)

  grand <- colMeans(scores)

  total_ss <- sum(
    rowSums(
      sweep(
        scores,
        2,
        grand,
        "-"
      )^2
    )
  )

  within_ss <- 0

  for (g in levels(cl)) {

    Xi <- scores[
      cl == g,
      ,
      drop = FALSE
    ]

    ctr <- colMeans(Xi)

    within_ss <- within_ss +
      sum(
        rowSums(
          sweep(
            Xi,
            2,
            ctr,
            "-"
          )^2
        )
      )
  }

  between_ss <- total_ss - within_ss

  (between_ss / (k - 1)) /
    (within_ss / (n - k))
}


davies_bouldin <- function(scores, cl) {

  scores <- as.matrix(scores)
  cl <- factor(cl)
  groups <- levels(cl)
  k <- length(groups)

  centroids <- do.call(
    rbind,
    lapply(
      groups,
      function(g) {
        colMeans(
          scores[
            cl == g,
            ,
            drop = FALSE
          ]
        )
      }
    )
  )

  scatter <- vapply(
    seq_along(groups),
    function(i) {

      Xi <- scores[
        cl == groups[i],
        ,
        drop = FALSE
      ]

      ctr <- centroids[i, ]

      mean(
        sqrt(
          rowSums(
            sweep(
              Xi,
              2,
              ctr,
              "-"
            )^2
          )
        )
      )
    },
    numeric(1)
  )

  centroid_dist <- as.matrix(
    stats::dist(
      centroids
    )
  )

  R <- matrix(
    NA_real_,
    nrow = k,
    ncol = k
  )

  for (i in seq_len(k)) {
    for (j in seq_len(k)) {

      if (i == j) {
        R[i, j] <- -Inf
      } else if (centroid_dist[i, j] == 0) {
        R[i, j] <- Inf
      } else {
        R[i, j] <-
          (scatter[i] + scatter[j]) /
          centroid_dist[i, j]
      }
    }
  }

  mean(
    apply(
      R,
      1,
      max
    )
  )
}

read_dpt_spectrum <- function(path) {

  # The DPT files supplied with this repository are plain-text
  # comma-separated exports with two numeric columns and no header:
  #
  #   wavenumber,absorbance
  #
  # Example:
  #   3655.71906,0.00208
  z <- utils::read.csv(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c(
      "wavenumber",
      "absorbance"
    )
  )

  if (ncol(z) != 2L) {
    stop(
      "Expected exactly two comma-separated columns in ",
      path
    )
  }

  out <- tibble::tibble(
    wavenumber = suppressWarnings(
      as.numeric(z$wavenumber)
    ),
    absorbance = suppressWarnings(
      as.numeric(z$absorbance)
    )
  ) |>
    dplyr::filter(
      is.finite(wavenumber),
      is.finite(absorbance)
    ) |>
    dplyr::arrange(
      dplyr::desc(wavenumber)
    )

  if (nrow(out) < 100L) {
    stop(
      "Too few spectral points parsed from ",
      path,
      ": ",
      nrow(out)
    )
  }

  out
}

# ============================================================
# PART A. ARCHIVED SINGLE CLEAN STRIP
# ============================================================

supp_ids <- c("etalon", "G", "N-1", "N-2", "N")

extract_supp <- function(id) {
  z <- spectra_long |>
    dplyr::filter(sample_id == id, retain_analysis_windows(wavenumber)) |>
    dplyr::arrange(dplyr::desc(wavenumber))
  if (nrow(z) != length(wn)) stop("Unexpected retained grid for ", id)
  if (!isTRUE(all.equal(as.numeric(z$wavenumber), wn, tolerance = 1e-12))) {
    stop("Wavenumber-grid mismatch for ", id)
  }
  as.numeric(z$absorbance)
}

X_sup_raw <- do.call(rbind, lapply(supp_ids, extract_supp))
rownames(X_sup_raw) <- supp_ids
colnames(X_sup_raw) <- colnames(X_snv)
X_sup_snv <- snv_matrix(X_sup_raw)

# Supplementary projection. Five PCs are used only for the
# reference-distance analysis; supplementary spectra do not define axes.
X_with_sup <- rbind(X_snv, X_sup_snv)
ind_sup <- (nrow(X_snv) + 1L):nrow(X_with_sup)
pca_ref <- FactoMineR::PCA(
  X_with_sup, scale.unit = FALSE, ncp = N_PC_REFERENCE,
  ind.sup = ind_sup, graph = FALSE
)
active5 <- as.matrix(pca_ref$ind$coord[, 1:N_PC_REFERENCE, drop = FALSE])
sup5 <- as.matrix(pca_ref$ind.sup$coord[, 1:N_PC_REFERENCE, drop = FALSE])
rownames(active5) <- ids
rownames(sup5) <- supp_ids

healthy_ids <- c("N-1", "N-2", "N")
healthy_centroid <- colMeans(sup5[healthy_ids, , drop = FALSE])
etalon_coord <- as.numeric(sup5["etalon", ])

dist_participant <- tibble::tibble(
  participant_id = ids,
  cluster = cluster_baseline,
  distance_to_healthy_centroid = sqrt(rowSums(sweep(active5, 2, healthy_centroid, "-")^2)),
  distance_to_archived_clean_strip = sqrt(rowSums(sweep(active5, 2, etalon_coord, "-")^2))
)

dist_cluster <- dist_participant |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n = dplyr::n(),
    healthy_mean = mean(distance_to_healthy_centroid),
    healthy_sd = stats::sd(distance_to_healthy_centroid),
    clean_strip_mean = mean(distance_to_archived_clean_strip),
    clean_strip_sd = stats::sd(distance_to_archived_clean_strip),
    .groups = "drop"
  )

archived_strip_snv <- as.numeric(X_sup_snv["etalon", ])
strip_similarity <- tibble::tibble(
  participant_id = ids,
  cluster = cluster_baseline,
  spearman_rho = apply(X_snv, 1, stats::cor, y = archived_strip_snv, method = "spearman")
)
strip_similarity_cluster <- strip_similarity |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n = dplyr::n(), mean_rho = mean(spearman_rho),
    median_rho = stats::median(spearman_rho), sd_rho = stats::sd(spearman_rho),
    min_rho = min(spearman_rho), max_rho = max(spearman_rho), .groups = "drop"
  )
strip_similarity_kw <- stats::kruskal.test(spearman_rho ~ factor(cluster), data = strip_similarity)

# Remove the component aligned with the archived clean-strip spectrum.
e <- archived_strip_snv - mean(archived_strip_snv)
X_centered <- t(apply(X_snv, 1, function(x) x - mean(x)))
beta <- apply(X_centered, 1, function(x) sum(x * e) / sum(e * e))
X_strip_resid <- X_centered - beta %*% t(e)
rownames(X_strip_resid) <- ids
colnames(X_strip_resid) <- colnames(X_snv)

pca_strip_resid <- FactoMineR::PCA(X_strip_resid, scale.unit = FALSE, ncp = 3, graph = FALSE)
scores_strip_resid <- as.matrix(pca_strip_resid$ind$coord[, 1:3, drop = FALSE])
cluster_strip_resid <- ward_k3(scores_strip_resid)
ARI_strip_resid <- mclust::adjustedRandIndex(cluster_baseline, cluster_strip_resid)

# Exclude the D1-dominant technical-overlap interval.
keep_no_D1 <- !(wn >= 870 & wn <= 1121)
X_no_D1 <- X_snv[, keep_no_D1, drop = FALSE]
pca_no_D1 <- FactoMineR::PCA(X_no_D1, scale.unit = FALSE, ncp = 3, graph = FALSE)
scores_no_D1 <- as.matrix(pca_no_D1$ind$coord[, 1:3, drop = FALSE])
cluster_no_D1 <- ward_k3(scores_no_D1)
ARI_no_D1 <- mclust::adjustedRandIndex(cluster_baseline, cluster_no_D1)

single_strip_summary <- tibble::tibble(
  metric = c(
    "Patient-strip Spearman KW p",
    "Strip-aligned component removal ARI vs baseline",
    "Strip-residualized cluster 1 n", "Strip-residualized cluster 2 n", "Strip-residualized cluster 3 n",
    "Excluding 870-1121 cm^-1 ARI vs baseline",
    "No-D1 cluster 1 n", "No-D1 cluster 2 n", "No-D1 cluster 3 n"
  ),
  value = c(
    strip_similarity_kw$p.value,
    ARI_strip_resid,
    as.integer(table(factor(cluster_strip_resid, levels = 1:3))),
    ARI_no_D1,
    as.integer(table(factor(cluster_no_D1, levels = 1:3)))
  )
)

# ============================================================
# PART B. FOUR DIFFERENT UNUSED STRIPS
# ============================================================

strip_files <- file.path(technical_dir, sprintf("etalon_%02d.dpt", 1:4))
if (!all(file.exists(strip_files))) {
  stop("Expected: ", paste(strip_files, collapse = ", "))
}

strip_list <- lapply(strip_files, read_dpt_spectrum)
names(strip_list) <- sprintf("strip_%02d", 1:4)
strip_list <- lapply(strip_list, function(z) {
  z |> dplyr::filter(retain_analysis_windows(wavenumber)) |>
    dplyr::arrange(dplyr::desc(wavenumber))
})

lens <- vapply(strip_list, nrow, integer(1))
if (length(unique(lens)) != 1L) stop("Four strips have different retained grid lengths.")
strip_wn <- strip_list[[1]]$wavenumber
for (i in seq_along(strip_list)) {
  if (!isTRUE(all.equal(strip_list[[i]]$wavenumber, strip_wn, tolerance = 1e-12))) {
    stop("Four strips have different wavenumber grids; no interpolation is performed.")
  }
}

X_strips <- do.call(rbind, lapply(strip_list, function(z) z$absorbance))
rownames(X_strips) <- names(strip_list)

pear <- stats::cor(t(X_strips), method = "pearson")
spear <- stats::cor(t(X_strips), method = "spearman")
pairs <- which(upper.tri(pear), arr.ind = TRUE)
pairwise_cor <- tibble::tibble(
  strip_1 = rownames(pear)[pairs[, 1]], strip_2 = colnames(pear)[pairs[, 2]],
  pearson_r = pear[pairs], spearman_rho = spear[pairs]
)

integrated <- tibble::tibble(
  strip = rownames(X_strips),
  integrated_absorbance = apply(X_strips, 1, function(y) integrated_two_windows(strip_wn, y))
)
integrated_cv <- cv_percent(integrated$integrated_absorbance)

pt_mean <- colMeans(X_strips)
pt_sd <- apply(X_strips, 2, stats::sd)
pt_cv <- 100 * pt_sd / abs(pt_mean)
pt_cv[!is.finite(pt_cv)] <- NA_real_
pointwise <- tibble::tibble(wavenumber = strip_wn, mean_absorbance = pt_mean, sd_absorbance = pt_sd, cv_percent = pt_cv)
pointwise_summary <- tibble::tibble(
  metric = c("Median pointwise CV (%)", "Mean pointwise CV (%)", "95th percentile pointwise CV (%)", "Maximum pointwise CV (%)"),
  value = c(stats::median(pt_cv, na.rm = TRUE), mean(pt_cv, na.rm = TRUE),
            stats::quantile(pt_cv, .95, na.rm = TRUE, names = FALSE), max(pt_cv, na.rm = TRUE))
)

domains <- tibble::tribble(
  ~domain, ~lower, ~upper,
  "D1", 876, 1127,
  "D2", 1150, 1313,
  "D3", 1382, 1413,
  "D4", 1477, 1706,
  "D5", 2942, 3092
)

domain_pairs <- dplyr::bind_rows(lapply(seq_len(nrow(domains)), function(i) {
  d <- domains[i, ]
  idx <- strip_wn >= d$lower & strip_wn <= d$upper
  p <- stats::cor(t(X_strips[, idx, drop = FALSE]), method = "pearson")
  s <- stats::cor(t(X_strips[, idx, drop = FALSE]), method = "spearman")
  ij <- which(upper.tri(p), arr.ind = TRUE)
  tibble::tibble(domain = d$domain, strip_1 = rownames(p)[ij[,1]], strip_2 = colnames(p)[ij[,2]],
                 pearson_r = p[ij], spearman_rho = s[ij])
}))

domain_summary <- domain_pairs |>
  dplyr::group_by(domain) |>
  dplyr::summarise(
    min_pearson = min(pearson_r), median_pearson = stats::median(pearson_r),
    min_spearman = min(spearman_rho), median_spearman = stats::median(spearman_rho),
    .groups = "drop"
  )

# Compare the archived reference with the four-strip set on the
# same instrumental grid. DPT text export stores wavenumbers to
# five decimal places, whereas the OPUS-derived pediatric grid
# retains additional decimal precision. Only this rounding-level
# difference is allowed; no interpolation is performed.
GRID_TOLERANCE_CM1 <- 1e-4

if (length(strip_wn) != length(wn)) {
  stop(
    "DPT and archived pediatric grids have different lengths; ",
    "no interpolation is performed."
  )
}

max_grid_difference <- max(
  abs(
    strip_wn - wn
  )
)

if (
  !is.finite(max_grid_difference) ||
  max_grid_difference > GRID_TOLERANCE_CM1
) {
  stop(
    "DPT grid differs from archived pediatric grid by more than ",
    GRID_TOLERANCE_CM1,
    " cm^-1. Maximum difference = ",
    signif(max_grid_difference, 6),
    " cm^-1; no interpolation is performed."
  )
}
archived_strip_raw <- as.numeric(X_sup_raw["etalon", ])
archived_vs_new <- tibble::tibble(
  strip = rownames(X_strips),
  pearson_r = apply(X_strips, 1, stats::cor, y = archived_strip_raw, method = "pearson"),
  spearman_rho = apply(X_strips, 1, stats::cor, y = archived_strip_raw, method = "spearman")
)
new_mean <- colMeans(X_strips)
archived_vs_mean <- tibble::tibble(
  pearson_r = stats::cor(archived_strip_raw, new_mean, method = "pearson"),
  spearman_rho = stats::cor(archived_strip_raw, new_mean, method = "spearman")
)

technical_summary <- tibble::tibble(
  metric = c(
    "Minimum pairwise Pearson", "Median pairwise Pearson",
    "Minimum pairwise Spearman", "Median pairwise Spearman",
    "Integrated absorbance CV (%)", "Median pointwise CV (%)",
    "95th percentile pointwise CV (%)",
    "Archived strip vs four-strip mean Pearson",
    "Archived strip vs four-strip mean Spearman"
  ),
  value = c(
    min(pairwise_cor$pearson_r), stats::median(pairwise_cor$pearson_r),
    min(pairwise_cor$spearman_rho), stats::median(pairwise_cor$spearman_rho),
    integrated_cv,
    pointwise_summary$value[pointwise_summary$metric == "Median pointwise CV (%)"],
    pointwise_summary$value[pointwise_summary$metric == "95th percentile pointwise CV (%)"],
    archived_vs_mean$pearson_r, archived_vs_mean$spearman_rho
  )
)


# ============================================================
# PART C. REFERENCE-INVARIANCE AND RESIDUAL-STRUCTURE TESTS
# ============================================================
#
# These analyses address two distinct technical questions:
#
#   1. Does the original "Cluster 1 is closest to clean strip"
#      observation depend on one particular archived strip?
#
#   2. Does removal of a substrate-aligned spectral direction
#      eliminate the multivariate organization of the pediatric
#      spectra, or does structured variation remain?
#
# They do NOT demonstrate that substrate contribution is absent.
# Instead, they distinguish:
#   - variability among the physical unused strips themselves;
#   - variable expression of a highly reproducible substrate-like
#     spectral component within patient spectra.
# ============================================================


# ------------------------------------------------------------
# 16. Put the four new strips on the pediatric analytical grid
# ------------------------------------------------------------
#
# The grid has already been shown to differ only by DPT text-export
# rounding (<= GRID_TOLERANCE_CM1), so columns can be assigned to the
# corresponding archived grid without interpolation.

colnames(X_strips) <- colnames(X_snv)

X_strips_snv <- snv_matrix(
  X_strips
)

consensus_four_strip_snv <- colMeans(
  X_strips_snv
)


# ------------------------------------------------------------
# 17. Project archived + four new strips into pediatric PCA space
# ------------------------------------------------------------
#
# All technical references are supplementary observations. They do
# not contribute to the PCA axes. Five PCs are retained here only to
# match the original descriptive reference-distance analysis.

X_projection_all <- rbind(
  X_snv,
  X_sup_snv,
  X_strips_snv
)

projection_sup_ids <- c(
  supp_ids,
  rownames(X_strips_snv)
)

ind_sup_all <- (
  nrow(X_snv) + 1L
):nrow(X_projection_all)

pca_all_strip_projection <- FactoMineR::PCA(
  X_projection_all,
  scale.unit = FALSE,
  ncp = N_PC_REFERENCE,
  ind.sup = ind_sup_all,
  graph = FALSE
)

active_projection_5 <- as.matrix(
  pca_all_strip_projection$ind$coord[
    ,
    1:N_PC_REFERENCE,
    drop = FALSE
  ]
)

supp_projection_5 <- as.matrix(
  pca_all_strip_projection$ind.sup$coord[
    ,
    1:N_PC_REFERENCE,
    drop = FALSE
  ]
)

rownames(active_projection_5) <- ids
rownames(supp_projection_5) <- projection_sup_ids

clean_strip_ids <- c(
  "etalon",
  rownames(X_strips_snv)
)

new_strip_ids <- rownames(
  X_strips_snv
)

new_strip_centroid_5 <- colMeans(
  supp_projection_5[
    new_strip_ids,
    ,
    drop = FALSE
  ]
)

clean_strip_projection_scores <- tibble::tibble(
  reference = clean_strip_ids,
  PC1 = supp_projection_5[
    clean_strip_ids,
    1
  ],
  PC2 = supp_projection_5[
    clean_strip_ids,
    2
  ],
  PC3 = supp_projection_5[
    clean_strip_ids,
    3
  ],
  PC4 = supp_projection_5[
    clean_strip_ids,
    4
  ],
  PC5 = supp_projection_5[
    clean_strip_ids,
    5
  ],
  distance_to_four_strip_centroid =
    sqrt(
      rowSums(
        sweep(
          supp_projection_5[
            clean_strip_ids,
            ,
            drop = FALSE
          ],
          2,
          new_strip_centroid_5,
          "-"
        )^2
      )
    )
)


# Baseline-cluster centroids in the same five-PC score space.

cluster_centroids_5 <- do.call(
  rbind,
  lapply(
    1:3,
    function(k) {
      colMeans(
        active_projection_5[
          cluster_baseline == k,
          ,
          drop = FALSE
        ]
      )
    }
  )
)

rownames(cluster_centroids_5) <- paste0(
  "Cluster ",
  1:3
)


# Distances from each baseline cluster centroid to:
#   - archived etalon,
#   - each of the four new strips,
#   - centroid of the four new strips.

reference_coords_list <- c(
  lapply(
    clean_strip_ids,
    function(id) {
      as.numeric(
        supp_projection_5[
          id,
          ,
          drop = TRUE
        ]
      )
    }
  ),
  list(
    new_four_strip_centroid =
      as.numeric(
        new_strip_centroid_5
      )
  )
)

names(reference_coords_list) <- c(
  clean_strip_ids,
  "new_four_strip_centroid"
)

strip_reference_cluster_distances <-
  dplyr::bind_rows(
    lapply(
      names(reference_coords_list),
      function(ref_name) {

        ref_coord <-
          reference_coords_list[[ref_name]]

        d <- sqrt(
          rowSums(
            sweep(
              cluster_centroids_5,
              2,
              ref_coord,
              "-"
            )^2
          )
        )

        tibble::tibble(
          reference = ref_name,
          cluster = 1:3,
          distance = as.numeric(d)
        )
      }
    )
  ) |>
  dplyr::group_by(
    reference
  ) |>
  dplyr::mutate(
    nearest_cluster =
      cluster[
        which.min(distance)
      ]
  ) |>
  dplyr::ungroup()


# ------------------------------------------------------------
# 18. Residualize against every clean-strip reference
# ------------------------------------------------------------

reference_vectors_snv <- rbind(
  archived_etalon =
    archived_strip_snv,
  X_strips_snv,
  mean_four_strips =
    consensus_four_strip_snv
)

residualization_objects <- lapply(
  seq_len(
    nrow(
      reference_vectors_snv
    )
  ),
  function(i) {

    ref_name <-
      rownames(
        reference_vectors_snv
      )[i]

    X_res <- remove_reference_component(
      X_snv,
      reference_vectors_snv[
        i,
        ,
        drop = TRUE
      ]
    )

    pca_res <- FactoMineR::PCA(
      X_res,
      scale.unit = FALSE,
      ncp = N_PC_BASELINE,
      graph = FALSE
    )

    scores_res <- as.matrix(
      pca_res$ind$coord[
        ,
        1:N_PC_BASELINE,
        drop = FALSE
      ]
    )

    cl_res <- ward_k3(
      scores_res
    )

    list(
      reference = ref_name,
      X = X_res,
      pca = pca_res,
      scores = scores_res,
      cluster = cl_res,
      ARI_vs_baseline =
        mclust::adjustedRandIndex(
          cluster_baseline,
          cl_res
        ),
      cumulative_variance_3PC =
        sum(
          pca_res$eig[
            1:N_PC_BASELINE,
            2
          ]
        )
    )
  }
)

names(residualization_objects) <-
  rownames(
    reference_vectors_snv
  )

archived_resid_cluster <-
  residualization_objects[["archived_etalon"]]$cluster


reference_residualization_summary <-
  dplyr::bind_rows(
    lapply(
      residualization_objects,
      function(obj) {

        tibble::tibble(
          reference = obj$reference,
          ARI_vs_baseline =
            obj$ARI_vs_baseline,
          ARI_vs_archived_residualized =
            mclust::adjustedRandIndex(
              archived_resid_cluster,
              obj$cluster
            ),
          cumulative_variance_3PC_percent =
            obj$cumulative_variance_3PC,
          cluster_sizes =
            paste(
              sort(
                as.integer(
                  table(
                    obj$cluster
                  )
                ),
                decreasing = TRUE
              ),
              collapse = "/"
            )
        )
      }
    )
  )


# Pairwise ARI among all residualized partitions gives a
# label-invariant measure of whether the conclusion depends on the
# particular unused-strip reference chosen.

resid_names <- names(
  residualization_objects
)

reference_residualization_pairwise_ARI <-
  dplyr::bind_rows(
    lapply(
      seq_along(resid_names),
      function(i) {

        if (i == length(resid_names)) {
          return(NULL)
        }

        dplyr::bind_rows(
          lapply(
            (i + 1):length(resid_names),
            function(j) {

              tibble::tibble(
                reference_1 =
                  resid_names[i],
                reference_2 =
                  resid_names[j],
                ARI =
                  mclust::adjustedRandIndex(
                    residualization_objects[[resid_names[i]]]$cluster,
                    residualization_objects[[resid_names[j]]]$cluster
                  )
              )
            }
          )
        )
      }
    )
  )


# ------------------------------------------------------------
# 19. k = 2-6 structure after consensus four-strip removal
# ------------------------------------------------------------
#
# This analysis does not "prove" a natural number of clusters.
# It asks whether removing the consensus substrate direction leaves
# nontrivial organization in the same parsimonious three-PC space.

mean_four_obj <-
  residualization_objects[["mean_four_strips"]]

mean_four_scores_3 <-
  mean_four_obj$scores

mean_four_dist_3 <-
  stats::dist(
    mean_four_scores_3,
    method = "euclidean"
  )

mean_four_hc <-
  stats::hclust(
    mean_four_dist_3,
    method = "ward.D2"
  )

mean_four_k_results <-
  dplyr::bind_rows(
    lapply(
      2:6,
      function(k) {

        cl_k <- stats::cutree(
          mean_four_hc,
          k = k
        )

        sil_k <- cluster::silhouette(
          cl_k,
          mean_four_dist_3
        )

        sizes_k <- as.integer(
          table(
            cl_k
          )
        )

        tibble::tibble(
          k = k,
          average_silhouette =
            mean(
              sil_k[
                ,
                "sil_width"
              ]
            ),
          calinski_harabasz =
            calinski_harabasz(
              mean_four_scores_3,
              cl_k
            ),
          davies_bouldin =
            davies_bouldin(
              mean_four_scores_3,
              cl_k
            ),
          smallest_cluster =
            min(
              sizes_k
            ),
          largest_cluster =
            max(
              sizes_k
            ),
          cluster_sizes =
            paste(
              sizes_k,
              collapse = "/"
            )
        )
      }
    )
  )

mean_four_k3_cluster <-
  stats::cutree(
    mean_four_hc,
    k = 3
  )

mean_four_k3_sil <- cluster::silhouette(
  mean_four_k3_cluster,
  mean_four_dist_3
)

mean_four_k3_diagnostics <- tibble::tibble(
  metric = c(
    "ARI vs baseline",
    "Average silhouette at k=3",
    "Negative silhouette count at k=3",
    "Three-PC cumulative variance (%)"
  ),
  value = c(
    mclust::adjustedRandIndex(
      cluster_baseline,
      mean_four_k3_cluster
    ),
    mean(
      mean_four_k3_sil[
        ,
        "sil_width"
      ]
    ),
    sum(
      mean_four_k3_sil[
        ,
        "sil_width"
      ] < 0
    ),
    mean_four_obj$
      cumulative_variance_3PC
  )
)


# Analytical figures; final manuscript assembly is deferred.
overlay <- tibble::tibble(
  wavenumber = rep(strip_wn, times = nrow(X_strips)),
  absorbance = as.numeric(t(X_strips)),
  strip = rep(rownames(X_strips), each = length(strip_wn))
)
p_overlay <- ggplot2::ggplot(overlay, ggplot2::aes(wavenumber, absorbance, group = strip)) +
  ggplot2::geom_line(linewidth = .55, alpha = .8) + ggplot2::scale_x_reverse() +
  ggplot2::labs(title = "Four independent unused-strip spectra", x = expression("Wavenumber (cm"^{-1}*")"), y = "Absorbance") +
  ggplot2::theme_minimal(base_size = 12)
ggplot2::ggsave(file.path(figures_dir, "technical_controls_four_strip_overlay.png"), p_overlay, width = 10, height = 5.5, dpi = 300)

p_cv <- ggplot2::ggplot(pointwise, ggplot2::aes(wavenumber, cv_percent)) +
  ggplot2::geom_line(linewidth = .6) + ggplot2::scale_x_reverse() +
  ggplot2::labs(title = "Pointwise technical variability across four unused strips", x = expression("Wavenumber (cm"^{-1}*")"), y = "Pointwise CV (%)") +
  ggplot2::theme_minimal(base_size = 12)
ggplot2::ggsave(file.path(figures_dir, "technical_controls_pointwise_CV.png"), p_cv, width = 10, height = 5.5, dpi = 300)

# Save numerical outputs.
readr::write_csv(dist_participant, file.path(tables_dir, "strip_archived_reference_distances_by_participant.csv"))
readr::write_csv(dist_cluster, file.path(tables_dir, "strip_archived_reference_distances_by_cluster.csv"))
readr::write_csv(strip_similarity, file.path(tables_dir, "strip_archived_patient_spearman_similarity.csv"))
readr::write_csv(strip_similarity_cluster, file.path(tables_dir, "strip_archived_similarity_by_cluster.csv"))
readr::write_csv(single_strip_summary, file.path(tables_dir, "strip_archived_sensitivity_summary.csv"))
readr::write_csv(pairwise_cor, file.path(tables_dir, "strip_four_controls_pairwise_correlations.csv"))
readr::write_csv(integrated, file.path(tables_dir, "strip_four_controls_integrated_absorbance.csv"))
readr::write_csv(pointwise, file.path(tables_dir, "strip_four_controls_pointwise_variability.csv"))
readr::write_csv(pointwise_summary, file.path(tables_dir, "strip_four_controls_pointwise_CV_summary.csv"))
readr::write_csv(domain_pairs, file.path(tables_dir, "strip_four_controls_domain_pairwise_correlations.csv"))
readr::write_csv(domain_summary, file.path(tables_dir, "strip_four_controls_domain_reproducibility_summary.csv"))
readr::write_csv(archived_vs_new, file.path(tables_dir, "strip_archived_vs_four_controls.csv"))
readr::write_csv(technical_summary, file.path(tables_dir, "strip_technical_reproducibility_summary.csv"))
readr::write_csv(clean_strip_projection_scores, file.path(tables_dir, "strip_all_clean_references_PCA_projection.csv"))
readr::write_csv(strip_reference_cluster_distances, file.path(tables_dir, "strip_reference_to_cluster_centroid_distances.csv"))
readr::write_csv(reference_residualization_summary, file.path(tables_dir, "strip_reference_residualization_sensitivity.csv"))
readr::write_csv(reference_residualization_pairwise_ARI, file.path(tables_dir, "strip_reference_residualization_pairwise_ARI.csv"))
readr::write_csv(mean_four_k_results, file.path(tables_dir, "strip_mean_four_residualized_k2_k6.csv"))
readr::write_csv(mean_four_k3_diagnostics, file.path(tables_dir, "strip_mean_four_residualized_k3_diagnostics.csv"))

saveRDS(
  list(
    settings = list(baseline_n_pc = 3L, reference_distance_n_pc = 5L, baseline_k = 3L, D1_exclusion = c(870,1121)),
    dist_cluster = dist_cluster,
    strip_similarity_cluster = strip_similarity_cluster,
    strip_similarity_kw = strip_similarity_kw,
    X_strip_resid = X_strip_resid,
    pca_strip_resid = pca_strip_resid,
    cluster_strip_resid = cluster_strip_resid,
    ARI_strip_resid = ARI_strip_resid,
    X_no_D1 = X_no_D1,
    pca_no_D1 = pca_no_D1,
    cluster_no_D1 = cluster_no_D1,
    ARI_no_D1 = ARI_no_D1,
    X_strips = X_strips,
    pairwise_cor = pairwise_cor,
    integrated = integrated,
    pointwise = pointwise,
    domain_summary = domain_summary,
    archived_vs_new = archived_vs_new,
    technical_summary = technical_summary,
    X_strips_snv = X_strips_snv,
    consensus_four_strip_snv = consensus_four_strip_snv,
    pca_all_strip_projection = pca_all_strip_projection,
    clean_strip_projection_scores = clean_strip_projection_scores,
    strip_reference_cluster_distances = strip_reference_cluster_distances,
    residualization_objects = residualization_objects,
    reference_residualization_summary = reference_residualization_summary,
    reference_residualization_pairwise_ARI = reference_residualization_pairwise_ARI,
    mean_four_k_results = mean_four_k_results,
    mean_four_k3_diagnostics = mean_four_k3_diagnostics
  ),
  file.path(processed_dir, "strip_reproducibility_objects.rds")
)

capture.output(sessionInfo(), file = file.path(tables_dir, "sessionInfo_strip_reproducibility.txt"))

cat("\n========================================\n")
cat("STRIP REPRODUCIBILITY ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("\nArchived reference distances (first five PCs):\n")
print(dist_cluster, n = Inf, width = Inf)
cat("\nArchived-strip similarity by cluster:\n")
print(strip_similarity_cluster, n = Inf, width = Inf)
cat("Archived-strip similarity KW p: ", signif(strip_similarity_kw$p.value, 5), "\n", sep = "")
cat("Strip-aligned residualization ARI: ", round(ARI_strip_resid, 3), "\n", sep = "")
cat("Residualized cluster sizes: ", paste(as.integer(table(factor(cluster_strip_resid, levels=1:3))), collapse="/"), "\n", sep="")
cat("Excluding 870-1121 cm^-1 ARI: ", round(ARI_no_D1, 3), "\n", sep = "")
cat("No-D1 cluster sizes: ", paste(as.integer(table(factor(cluster_no_D1, levels=1:3))), collapse="/"), "\n", sep="")
cat("\nFour unused strips:\n")
cat(
  "Retained DPT points per strip: ",
  ncol(X_strips),
  "\n",
  sep = ""
)
cat(
  "Maximum DPT-vs-archived grid difference: ",
  signif(max_grid_difference, 6),
  " cm^-1\n",
  sep = ""
)
cat("Pairwise Pearson range: ", paste(signif(range(pairwise_cor$pearson_r), 8), collapse=" to "), "\n", sep="")
cat("Pairwise Spearman range: ", paste(signif(range(pairwise_cor$spearman_rho), 8), collapse=" to "), "\n", sep="")
cat("Integrated absorbance CV: ", round(integrated_cv, 3), "%\n", sep="")
cat("Median pointwise CV: ", round(stats::median(pt_cv, na.rm=TRUE), 3), "%\n", sep="")
cat("95th percentile pointwise CV: ", round(stats::quantile(pt_cv, .95, na.rm=TRUE), 3), "%\n", sep="")
cat("\nDomain-specific reproducibility:\n")
print(domain_summary, n = Inf, width = Inf)
cat("\nArchived strip vs each new strip:\n")
print(archived_vs_new, n = Inf, width = Inf)
cat("\nArchived strip vs four-strip mean:\n")
print(archived_vs_mean, n = Inf, width = Inf)

cat("\nAll clean-strip references projected into pediatric PCA space:\n")
print(
  clean_strip_projection_scores,
  n = Inf,
  width = Inf
)

cat("\nDistance from each clean-strip reference to baseline cluster centroids (5 PCs):\n")
print(
  strip_reference_cluster_distances,
  n = Inf,
  width = Inf
)

cat("\nResidualization using each clean-strip reference:\n")
print(
  reference_residualization_summary,
  n = Inf,
  width = Inf
)

cat("\nPairwise ARI among residualized partitions:\n")
print(
  reference_residualization_pairwise_ARI,
  n = Inf,
  width = Inf
)

cat("\nConsensus four-strip removal: k = 2-6 internal indices:\n")
print(
  mean_four_k_results,
  n = Inf,
  width = Inf
)

cat("\nConsensus four-strip removal: k = 3 diagnostics:\n")
print(
  mean_four_k3_diagnostics,
  n = Inf,
  width = Inf
)

cat("========================================\n")
