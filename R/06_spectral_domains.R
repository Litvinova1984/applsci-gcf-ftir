# ============================================================
# 06_spectral_domains.R
#
# Pediatric GCF FTIR: PLS-DA/VIP spectral-domain analysis
#
# This script:
#   1. fits a three-component PLS-DA model to the frozen pediatric
#      SNV-normalized spectra and baseline cluster labels;
#   2. summarizes VIP scores as the maximum VIP across components;
#   3. identifies contiguous VIP-enriched spectral regions (VIP > 1);
#   4. computes participant-level scores for the fixed D1-D5 domains;
#   5. performs Kruskal-Wallis and pairwise Wilcoxon-BH tests;
#   6. derives cluster-level standardized domain signatures;
#   7. prepares median raw absorbance spectra + IQR for later figures.
#
# PLS-DA is used as a secondary supervised characterization of
# spectral variables associated with the already-defined baseline
# clusters. It is not used to select or validate the cluster solution.
# ============================================================

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "00_utils.R"))

required_pkgs <- c(
  "mixOmics", "tibble", "dplyr", "tidyr", "readr", "rstatix"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

set.seed(SEED)

processed_dir <- file.path("data", "pediatric", "processed")
tables_dir <- "tables"

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

x_snv_file <- file.path(processed_dir, "pediatric_X_snv.rds")
x_raw_file <- file.path(processed_dir, "pediatric_X_active_absorbance.rds")
wn_file <- file.path(processed_dir, "pediatric_wavenumbers.rds")
cluster_file <- file.path(processed_dir, "pediatric_clusters_baseline.rds")

for (f in c(x_snv_file, x_raw_file, wn_file, cluster_file)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}

X_snv <- as.matrix(readRDS(x_snv_file))
X_raw <- as.matrix(readRDS(x_raw_file))
wn_precise <- as.numeric(readRDS(wn_file))
cluster_raw <- readRDS(cluster_file)

if (!identical(dim(X_snv), c(25L, 793L))) {
  stop(
    "Expected X_snv dimensions 25 x 793; obtained ",
    paste(dim(X_snv), collapse = " x ")
  )
}

if (!identical(dim(X_raw), dim(X_snv))) {
  stop("Raw and SNV matrices have different dimensions.")
}

if (length(wn_precise) != ncol(X_snv)) {
  stop("Wavenumber vector length does not match spectral matrices.")
}

if (is.null(rownames(X_snv)) || is.null(rownames(X_raw))) {
  stop("Spectral matrices must have participant IDs as row names.")
}

if (!setequal(rownames(X_snv), rownames(X_raw))) {
  stop("Raw and SNV matrices do not contain the same participant IDs.")
}

ids <- rownames(X_snv)
X_raw <- X_raw[ids, , drop = FALSE]

if (is.null(names(cluster_raw))) {
  if (length(cluster_raw) != length(ids)) {
    stop("Unnamed baseline cluster vector has wrong length.")
  }
  cluster_baseline <- as.integer(cluster_raw)
  names(cluster_baseline) <- ids
} else {
  if (!setequal(names(cluster_raw), ids)) {
    stop("Baseline cluster IDs do not match pediatric spectra.")
  }

  cluster_raw <- cluster_raw[ids]

  if (is.factor(cluster_raw)) {
    cluster_baseline <- as.integer(as.character(cluster_raw))
  } else {
    cluster_baseline <- as.integer(cluster_raw)
  }

  names(cluster_baseline) <- ids
}

cluster_factor <- factor(cluster_baseline, levels = 1:3)

if (!identical(as.integer(table(cluster_factor)), c(15L, 4L, 6L))) {
  stop("Expected baseline cluster sizes 15 / 4 / 6.")
}

# Processed matrices use two-decimal wavenumber labels.
# The precise OPUS-derived grid is retained separately for raw spectra.
if (is.null(colnames(X_snv))) {
  stop("X_snv must have wavenumber column names.")
}

wn_domain <- suppressWarnings(
  as.numeric(sub("^wn_", "", colnames(X_snv)))
)

if (anyNA(wn_domain)) {
  stop("Could not parse wavenumbers from X_snv column names.")
}

if (max(abs(wn_domain - wn_precise)) > 0.01) {
  stop("Processed and precise wavenumber grids differ unexpectedly.")
}

# ------------------------------------------------------------
# Three-component PLS-DA and VIP profile
# ------------------------------------------------------------

N_PLS_COMPONENTS <- 3L

plsda_fit <- mixOmics::plsda(
  X_snv,
  cluster_factor,
  ncomp = N_PLS_COMPONENTS
)

vip_matrix <- mixOmics::vip(plsda_fit)

if (nrow(vip_matrix) != ncol(X_snv)) {
  stop("VIP matrix does not match number of spectral variables.")
}

if (ncol(vip_matrix) < N_PLS_COMPONENTS) {
  stop("VIP matrix contains fewer than three components.")
}

vip_summary <- tibble::tibble(
  variable = rownames(vip_matrix),
  wavenumber = wn_domain,
  VIP_comp1 = as.numeric(vip_matrix[, 1]),
  VIP_comp2 = as.numeric(vip_matrix[, 2]),
  VIP_comp3 = as.numeric(vip_matrix[, 3]),
  VIP_max = apply(
    vip_matrix[, 1:N_PLS_COMPONENTS, drop = FALSE],
    1,
    max
  )
) |>
  dplyr::arrange(dplyr::desc(VIP_max))

# Variables with summary VIP > 1 are ordered by wavenumber.
# A gap > 5 cm^-1 starts a new region.
vip_gt1 <- vip_summary |>
  dplyr::filter(VIP_max > 1) |>
  dplyr::arrange(wavenumber) |>
  dplyr::mutate(
    gap_cm1 = c(0, diff(wavenumber)),
    vip_region = cumsum(gap_cm1 > 5)
  )

vip_regions <- vip_gt1 |>
  dplyr::group_by(vip_region) |>
  dplyr::summarise(
    start = min(wavenumber),
    end = max(wavenumber),
    n_variables = dplyr::n(),
    mean_VIP = mean(VIP_max),
    max_VIP = max(VIP_max),
    .groups = "drop"
  ) |>
  dplyr::arrange(start)

# ------------------------------------------------------------
# Fixed D1-D5 spectral domains
# ------------------------------------------------------------

domain_definitions <- tibble::tribble(
  ~domain, ~lower, ~upper,
  "D1",  876, 1127,
  "D2", 1150, 1313,
  "D3", 1382, 1413,
  "D4", 1477, 1706,
  "D5", 2942, 3092
)

domain_list <- setNames(
  lapply(
    seq_len(nrow(domain_definitions)),
    function(i) {
      c(domain_definitions$lower[i], domain_definitions$upper[i])
    }
  ),
  domain_definitions$domain
)

vip_domain_coverage <- domain_definitions |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_VIP_gt1 = sum(
      vip_gt1$wavenumber >= lower &
        vip_gt1$wavenumber <= upper
    ),
    max_VIP = {
      vals <- vip_gt1$VIP_max[
        vip_gt1$wavenumber >= lower &
          vip_gt1$wavenumber <= upper
      ]
      if (length(vals) == 0L) NA_real_ else max(vals)
    },
    mean_VIP = {
      vals <- vip_gt1$VIP_max[
        vip_gt1$wavenumber >= lower &
          vip_gt1$wavenumber <= upper
      ]
      if (length(vals) == 0L) NA_real_ else mean(vals)
    }
  ) |>
  dplyr::ungroup()

# ------------------------------------------------------------
# Participant-level domain scores
# ------------------------------------------------------------

domain_scores <- tibble::tibble(participant_id = ids)
domain_grid_summary <- vector("list", length(domain_list))

for (i in seq_along(domain_list)) {
  domain_name <- names(domain_list)[i]
  rng <- domain_list[[i]]

  idx <- which(
    wn_domain >= rng[1] &
      wn_domain <= rng[2]
  )

  if (length(idx) == 0L) {
    stop("No spectral variables found for ", domain_name, ".")
  }

  domain_scores[[domain_name]] <- rowMeans(
    X_snv[, idx, drop = FALSE]
  )

  domain_grid_summary[[i]] <- tibble::tibble(
    domain = domain_name,
    lower = rng[1],
    upper = rng[2],
    n_wavenumbers = length(idx),
    minimum_retained_wavenumber = min(wn_domain[idx]),
    maximum_retained_wavenumber = max(wn_domain[idx])
  )
}

domain_grid_summary <- dplyr::bind_rows(domain_grid_summary)
domain_scores$cluster <- cluster_baseline

domain_long <- domain_scores |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(names(domain_list)),
    names_to = "domain",
    values_to = "value"
  ) |>
  dplyr::mutate(
    domain = factor(domain, levels = names(domain_list)),
    cluster = factor(cluster, levels = 1:3)
  )

# ------------------------------------------------------------
# Domain-level statistical tests
# ------------------------------------------------------------

domain_kw <- domain_long |>
  dplyr::group_by(domain) |>
  rstatix::kruskal_test(value ~ cluster) |>
  dplyr::ungroup() |>
  dplyr::arrange(domain)

domain_pairwise <- domain_long |>
  dplyr::group_by(domain) |>
  rstatix::pairwise_wilcox_test(
    value ~ cluster,
    p.adjust.method = "BH"
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(domain, group1, group2)

# ------------------------------------------------------------
# Cluster-level standardized domain signatures
# ------------------------------------------------------------

cluster_domain_means <- domain_long |>
  dplyr::group_by(cluster, domain) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_value = mean(value, na.rm = TRUE),
    median_value = stats::median(value, na.rm = TRUE),
    sd_value = stats::sd(value, na.rm = TRUE),
    .groups = "drop"
  )

cluster_domain_signatures <- cluster_domain_means |>
  dplyr::group_by(domain) |>
  dplyr::mutate(
    standardized_effect = as.numeric(scale(mean_value)),
    direction = dplyr::if_else(
      standardized_effect >= 0,
      "higher",
      "lower"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(cluster, domain)

# ------------------------------------------------------------
# Median raw absorbance spectra + IQR by cluster
# ------------------------------------------------------------

raw_spectral_summary <- dplyr::bind_rows(
  lapply(
    1:3,
    function(k) {
      Xk <- X_raw[
        cluster_baseline == k,
        ,
        drop = FALSE
      ]

      tibble::tibble(
        cluster = k,
        wavenumber = wn_precise,
        median_absorbance = apply(
          Xk,
          2,
          stats::median,
          na.rm = TRUE
        ),
        q25_absorbance = apply(
          Xk,
          2,
          stats::quantile,
          probs = 0.25,
          na.rm = TRUE,
          names = FALSE
        ),
        q75_absorbance = apply(
          Xk,
          2,
          stats::quantile,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE
        ),
        mean_absorbance = colMeans(
          Xk,
          na.rm = TRUE
        ),
        n = nrow(Xk)
      )
    }
  )
)

domain_summary_for_reporting <- cluster_domain_means |>
  dplyr::left_join(
    domain_kw |>
      dplyr::select(domain, statistic, p),
    by = "domain"
  ) |>
  dplyr::left_join(
    domain_definitions,
    by = "domain"
  ) |>
  dplyr::arrange(domain, cluster)

# ------------------------------------------------------------
# Save tables and reusable objects
# ------------------------------------------------------------

readr::write_csv(
  vip_summary,
  file.path(tables_dir, "spectral_domains_VIP_all_variables.csv")
)

readr::write_csv(
  vip_gt1,
  file.path(tables_dir, "spectral_domains_VIP_gt1_variables.csv")
)

readr::write_csv(
  vip_regions,
  file.path(tables_dir, "spectral_domains_VIP_regions.csv")
)

readr::write_csv(
  vip_domain_coverage,
  file.path(tables_dir, "spectral_domains_VIP_domain_coverage.csv")
)

readr::write_csv(
  domain_grid_summary,
  file.path(tables_dir, "spectral_domains_grid_summary.csv")
)

readr::write_csv(
  domain_scores,
  file.path(tables_dir, "spectral_domains_participant_scores.csv")
)

readr::write_csv(
  domain_long,
  file.path(tables_dir, "spectral_domains_participant_scores_long.csv")
)

readr::write_csv(
  domain_kw,
  file.path(tables_dir, "spectral_domains_kruskal_wallis.csv")
)

readr::write_csv(
  domain_pairwise,
  file.path(tables_dir, "spectral_domains_pairwise_wilcoxon_BH.csv")
)

readr::write_csv(
  cluster_domain_means,
  file.path(tables_dir, "spectral_domains_cluster_means.csv")
)

readr::write_csv(
  cluster_domain_signatures,
  file.path(tables_dir, "spectral_domains_cluster_standardized_signatures.csv")
)

readr::write_csv(
  raw_spectral_summary,
  file.path(tables_dir, "spectral_domains_raw_absorbance_cluster_summary.csv")
)

readr::write_csv(
  domain_summary_for_reporting,
  file.path(tables_dir, "spectral_domains_reporting_summary.csv")
)

spectral_domain_objects <- list(
  settings = list(
    n_pls_components = N_PLS_COMPONENTS,
    domain_definitions = domain_definitions,
    baseline_cluster_sizes = as.integer(table(cluster_factor))
  ),
  plsda_fit = plsda_fit,
  vip_matrix = vip_matrix,
  vip_summary = vip_summary,
  vip_gt1 = vip_gt1,
  vip_regions = vip_regions,
  vip_domain_coverage = vip_domain_coverage,
  domain_scores = domain_scores,
  domain_long = domain_long,
  domain_kw = domain_kw,
  domain_pairwise = domain_pairwise,
  cluster_domain_means = cluster_domain_means,
  cluster_domain_signatures = cluster_domain_signatures,
  raw_spectral_summary = raw_spectral_summary
)

saveRDS(
  spectral_domain_objects,
  file.path(processed_dir, "spectral_domain_objects.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(tables_dir, "sessionInfo_spectral_domains.txt")
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n")
cat("========================================\n")
cat("SPECTRAL DOMAIN ANALYSIS COMPLETE\n")
cat("========================================\n")
cat(
  "Input matrix: ",
  nrow(X_snv),
  " x ",
  ncol(X_snv),
  "\n",
  sep = ""
)
cat(
  "Baseline cluster sizes: ",
  paste(as.integer(table(cluster_factor)), collapse = "/"),
  "\n",
  sep = ""
)
cat("PLS-DA components: ", N_PLS_COMPONENTS, "\n", sep = "")
cat("Variables with VIP > 1: ", nrow(vip_gt1), "\n", sep = "")

cat("\nVIP-enriched spectral regions:\n")
print(vip_regions, n = Inf, width = Inf)

cat("\nVIP > 1 coverage within D1-D5:\n")
print(vip_domain_coverage, n = Inf, width = Inf)

cat("\nKruskal-Wallis tests:\n")
print(domain_kw, n = Inf, width = Inf)

cat("\nPairwise Wilcoxon tests with BH correction:\n")
print(domain_pairwise, n = Inf, width = Inf)

cat("\nSigned standardized cluster-domain signatures:\n")
print(
  cluster_domain_signatures |>
    dplyr::select(
      cluster,
      domain,
      mean_value,
      standardized_effect,
      direction
    ),
  n = Inf,
  width = Inf
)

cat("========================================\n")
