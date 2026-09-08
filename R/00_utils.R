# ============================================================
# 00_utils.R
# Shared utilities and project constants
# ============================================================

SEED <- 4489628L

analysis_windows <- list(
  high = c(2800, 3400),
  low  = c(870, 1800)
)

spectral_domains <- data.frame(
  domain = c("D1", "D2", "D3", "D4", "D5"),
  xmin   = c(876, 1150, 1382, 1477, 2942),
  xmax   = c(1127, 1313, 1413, 1706, 3092)
)


# ------------------------------------------------------------
# Standard normal variate
# ------------------------------------------------------------

snv <- function(x) {
  
  x <- as.numeric(x)
  
  s <- stats::sd(
    x,
    na.rm = TRUE
  )
  
  if (!is.finite(s) || s == 0) {
    stop("SNV cannot be calculated: SD is zero or non-finite.")
  }
  
  (
    x -
      mean(x, na.rm = TRUE)
  ) / s
}


snv_matrix <- function(X) {
  
  X <- as.matrix(X)
  
  out <- t(
    apply(
      X,
      1,
      snv
    )
  )
  
  rownames(out) <- rownames(X)
  colnames(out) <- colnames(X)
  
  out
}


# ------------------------------------------------------------
# Historical spectral-column naming
# ------------------------------------------------------------

make_wn_names <- function(x) {
  
  x <- round(
    as.numeric(x),
    2
  )
  
  x_chr <- sprintf(
    "%.2f",
    x
  )
  
  x_chr <- sub(
    "0+$",
    "",
    x_chr
  )
  
  x_chr <- sub(
    "\\.$",
    "",
    x_chr
  )
  
  paste0(
    "wn_",
    x_chr
  )
}


extract_wavenumber <- function(x) {
  
  as.numeric(
    sub(
      "^wn_",
      "",
      x
    )
  )
}


# ------------------------------------------------------------
# Analysis-window filter
# ------------------------------------------------------------

retain_analysis_windows <- function(wavenumber) {
  
  (
    wavenumber >= analysis_windows$high[1] &
      wavenumber <= analysis_windows$high[2]
  ) |
    (
      wavenumber >= analysis_windows$low[1] &
        wavenumber <= analysis_windows$low[2]
    )
}


# ------------------------------------------------------------
# Duplicate check
# ------------------------------------------------------------

assert_no_duplicates <- function(x, label = "values") {
  
  if (anyDuplicated(x)) {
    stop(
      "Duplicate ",
      label,
      " detected."
    )
  }
  
  invisible(TRUE)
}