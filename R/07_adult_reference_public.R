# ============================================================
# 07_adult_reference_public.R
#
# Public reproducibility module for the archived adult GCF reference.
#
# Public inputs contain:
#   - a de-identified harmonized adult SNV spectral matrix (18 x 793);
#   - the retained wavenumber grid;
#   - frozen aggregate clinical-axis and historical-axis summaries.
#
# Individual-level adult clinical metadata, source labels, and broader-window
# source data are intentionally not distributed. This module independently
# recomputes the harmonized adult PCA and D1-D5 loading summaries and checks
# them against frozen aggregate reference outputs. It performs no new
# inferential analysis.
# ============================================================

if (!file.exists(file.path("R", "00_utils.R"))) {
  stop("Run this script from the project root.")
}
source(file.path("R", "00_utils.R"))

required_pkgs <- c("FactoMineR", "tibble", "dplyr", "tidyr", "readr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "))
}

set.seed(SEED)

processed_dir <- file.path("data", "adult", "processed")
aggregate_dir <- file.path("data", "adult", "aggregate")
tables_dir <- "tables"
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

x_file <- file.path(processed_dir, "adult_reference_X_snv.rds")
wn_file <- file.path(processed_dir, "adult_reference_wavenumbers.rds")

aggregate_files <- c(
  sign = "adult_PCA_axis_sign_convention.csv",
  clinical = "adult_PCA_selected_axis_clinical_summary.csv",
  category = "adult_PCA_selected_axis_category_means.csv",
  matches = "adult_PCA_selected_axis_matches.csv",
  correlations = "adult_PCA_historical_harmonized_score_correlations.csv",
  aligned = "adult_PCA_variable_coordinates_aligned_PC2_PC4_PC5.csv",
  domain_map = "adult_domain_clinical_mapping.csv"
)

required_files <- c(
  x_file,
  wn_file,
  file.path(aggregate_dir, unname(aggregate_files))
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0L) {
  stop("Required public adult input(s) missing:\n", paste(missing, collapse = "\n"))
}

read_agg <- function(key) {
  readr::read_csv(
    file.path(aggregate_dir, aggregate_files[[key]]),
    show_col_types = FALSE
  )
}

X <- as.matrix(readRDS(x_file))
wn <- as.numeric(readRDS(wn_file))

sign_tbl <- read_agg("sign")
clinical_summary <- read_agg("clinical")
category_means <- read_agg("category")
selected_matches <- read_agg("matches")
score_correlations <- read_agg("correlations")
frozen_aligned <- read_agg("aligned")
frozen_domain_map <- read_agg("domain_map")

if (!identical(dim(X), c(18L, 793L))) stop("Expected adult X_snv = 18 x 793.")
if (length(wn) != 793L) stop("Expected 793 retained adult wavenumbers.")
if (any(!is.finite(X)) || any(!is.finite(wn))) stop("Non-finite adult values detected.")
if (max(abs(rowMeans(X))) > 1e-10) stop("Adult SNV row-mean check failed.")
if (max(abs(apply(X, 1, stats::sd) - 1)) > 1e-10) stop("Adult SNV row-SD check failed.")

pca <- FactoMineR::PCA(
  X,
  scale.unit = FALSE,
  ncp = min(10L, nrow(X) - 1L, ncol(X)),
  graph = FALSE
)

adult_eigenvalues <- tibble::tibble(
  PC = seq_len(nrow(pca$eig)),
  eigenvalue = pca$eig[, 1],
  variance_percent = pca$eig[, 2],
  cumulative_variance_percent = pca$eig[, 3]
)

raw_coord <- tibble::as_tibble(
  as.data.frame(pca$var$coord),
  rownames = "variable"
) |>
  dplyr::mutate(wavenumber = wn, .after = variable) |>
  dplyr::transmute(
    wavenumber,
    variable,
    PC1 = Dim.1,
    PC2 = Dim.2,
    PC3 = Dim.3,
    PC4 = Dim.4,
    PC5 = Dim.5
  )

required_sign_pcs <- c("PC2", "PC4", "PC5")
if (!all(required_sign_pcs %in% sign_tbl$PC)) {
  stop("Frozen adult aggregate inputs lack PC2/PC4/PC5 sign conventions.")
}
flip <- setNames(sign_tbl$flip, sign_tbl$PC)

aligned_coord <- raw_coord |>
  dplyr::transmute(
    wavenumber,
    variable,
    PC2 = PC2 * flip[["PC2"]],
    PC4 = PC4 * flip[["PC4"]],
    PC5 = PC5 * flip[["PC5"]]
  )

# Exact-grid integrity check against the frozen aggregate loading export.
frozen_aligned_check <- frozen_aligned |>
  dplyr::arrange(dplyr::desc(wavenumber))
aligned_check <- aligned_coord |>
  dplyr::arrange(dplyr::desc(wavenumber))

if (nrow(frozen_aligned_check) != 793L) {
  stop("Frozen aligned adult loading reference must contain 793 rows.")
}
for (pc in required_sign_pcs) {
  r <- stats::cor(aligned_check[[pc]], frozen_aligned_check[[pc]], method = "pearson")
  if (!is.finite(r) || abs(r) < 0.999999) {
    stop("Recomputed adult ", pc, " loadings do not match the frozen aggregate reference.")
  }
}

domain_defs <- tibble::tribble(
  ~domain, ~lower, ~upper,
  "D1",  876, 1127,
  "D2", 1150, 1313,
  "D3", 1382, 1413,
  "D4", 1477, 1706,
  "D5", 2942, 3092
)

domain_axis_map <- tibble::tribble(
  ~domain, ~PC,
  "D1", "PC2",
  "D2", "PC2",
  "D3", "PC4",
  "D4", "PC4",
  "D5", "PC5"
)

loading_long <- aligned_coord |>
  tidyr::pivot_longer(
    cols = c(PC2, PC4, PC5),
    names_to = "PC",
    values_to = "loading"
  )

domain_all <- dplyr::bind_rows(
  lapply(seq_len(nrow(domain_defs)), function(i) {
    d <- domain_defs[i, ]
    loading_long |>
      dplyr::filter(wavenumber >= d$lower, wavenumber <= d$upper) |>
      dplyr::group_by(PC) |>
      dplyr::summarise(
        domain = d$domain,
        lower = d$lower,
        upper = d$upper,
        n_points = dplyr::n(),
        mean_loading = mean(loading),
        median_loading = stats::median(loading),
        min_loading = min(loading),
        max_loading = max(loading),
        max_abs_loading = loading[which.max(abs(loading))],
        prop_positive = mean(loading > 0),
        prop_negative = mean(loading < 0),
        .groups = "drop"
      )
  })
) |>
  dplyr::mutate(
    direction_class = dplyr::case_when(
      prop_positive >= 0.75 ~ "predominantly positive",
      prop_negative >= 0.75 ~ "predominantly negative",
      TRUE ~ "mixed"
    )
  ) |>
  dplyr::arrange(domain, PC)

domain_selected <- domain_axis_map |>
  dplyr::left_join(domain_all, by = c("domain", "PC"))

# Attach only frozen aggregate clinical context; no individual adult records.
domain_mapping <- domain_selected |>
  dplyr::left_join(
    frozen_domain_map |>
      dplyr::select(domain, PC, axis_context, domain_level_interpretation),
    by = c("domain", "PC")
  )

# Check selected-domain numerical summaries against the frozen aggregate map.
check_cols <- c(
  "n_points", "mean_loading", "median_loading",
  "prop_positive", "prop_negative"
)
for (i in seq_len(nrow(domain_mapping))) {
  d <- domain_mapping$domain[i]
  pc <- domain_mapping$PC[i]
  ref <- frozen_domain_map |>
    dplyr::filter(domain == d, PC == pc)
  if (nrow(ref) != 1L) stop("Frozen adult domain map has no unique row for ", d, "/", pc)
  for (nm in check_cols) {
    if (!isTRUE(all.equal(
      as.numeric(domain_mapping[[nm]][i]),
      as.numeric(ref[[nm]][1]),
      tolerance = 1e-10
    ))) {
      stop("Adult domain check failed for ", d, "/", pc, ": ", nm)
    }
  }
  if (!identical(domain_mapping$direction_class[i], ref$direction_class[1])) {
    stop("Adult direction-class check failed for ", d, "/", pc)
  }
}

# Frozen selected-axis mapping must match the manuscript specification.
expected_matches <- tibble::tribble(
  ~historical_PC, ~harmonized_PC,
  "PC1", "PC2",
  "PC2", "PC1",
  "PC4", "PC4",
  "PC5", "PC5"
)
observed_matches <- selected_matches |>
  dplyr::select(historical_PC, harmonized_PC) |>
  dplyr::arrange(historical_PC)
if (!identical(
  as.data.frame(observed_matches),
  as.data.frame(expected_matches |> dplyr::arrange(historical_PC))
)) {
  stop("Frozen historical-to-harmonized adult axis mapping is unexpected.")
}

# Publication/reviewer-facing aggregate outputs.
readr::write_csv(adult_eigenvalues, file.path(tables_dir, "adult_PCA_eigenvalues.csv"))
readr::write_csv(raw_coord, file.path(tables_dir, "adult_PCA_variable_coordinates_raw_PC1_PC5.csv"))
readr::write_csv(aligned_coord, file.path(tables_dir, "adult_PCA_variable_coordinates_aligned_PC2_PC4_PC5.csv"))
readr::write_csv(sign_tbl, file.path(tables_dir, "adult_PCA_axis_sign_convention.csv"))
readr::write_csv(domain_all, file.path(tables_dir, "adult_PCA_domain_loading_summary_all_selected_axes.csv"))
readr::write_csv(domain_selected, file.path(tables_dir, "adult_PCA_domain_loading_summary_selected_axes.csv"))
readr::write_csv(clinical_summary, file.path(tables_dir, "adult_PCA_selected_axis_clinical_summary.csv"))
readr::write_csv(category_means, file.path(tables_dir, "adult_PCA_selected_axis_category_means.csv"))
readr::write_csv(selected_matches, file.path(tables_dir, "adult_PCA_selected_axis_matches.csv"))
readr::write_csv(score_correlations, file.path(tables_dir, "adult_PCA_historical_harmonized_score_correlations.csv"))
readr::write_csv(domain_mapping, file.path(tables_dir, "adult_domain_clinical_mapping.csv"))

# Minimal public object consumed by script 08 and the Quarto report.
# It contains aggregate summaries only, never individual adult clinical data.
public_obj <- list(
  settings = list(
    n_reference_spectra = nrow(X),
    retained_variables = ncol(X),
    analytical_windows = c("3400-2800", "1800-870"),
    preprocessing = "row-wise SNV; PCA scale.unit = FALSE"
  ),
  axis_matching = list(
    selected_matches = selected_matches,
    score_correlations = score_correlations
  ),
  adult_clinical_eta2 = clinical_summary,
  adult_category_means = category_means,
  adult_domain_clinical_mapping = domain_mapping,
  adult_varcoord_aligned = aligned_coord
)

saveRDS(pca, file.path(processed_dir, "adult_reference_pca.rds"))
saveRDS(public_obj, file.path(processed_dir, "adult_reference_analysis_objects.rds"))

capture.output(sessionInfo(), file = file.path(tables_dir, "sessionInfo_adult_reference.txt"))

cat("\n========================================\n")
cat("PUBLIC ADULT REFERENCE CHECK COMPLETE\n")
cat("========================================\n")
cat("Reference spectra: 18\n")
cat("Retained variables: 793\n")
cat("Harmonized windows: 3400-2800 + 1800-870 cm^-1\n")
cat("Historical axis mapping (frozen aggregate): PC1->PC2; PC2->PC1; PC4->PC4; PC5->PC5\n")
cat("Domain map: D1/D2->PC2; D3/D4->PC4; D5->PC5\n")
cat("Individual-level adult clinical metadata: not used or released\n")
cat("========================================\n")
