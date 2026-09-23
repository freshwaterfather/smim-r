#' Noise level, arrival time and detection-limit tail cutoff for a breakthrough curve
#'
#' Rule-based helpers applied identically to every curve. They only *locate* times; the
#' raw data are what get fitted (no clamping, imputation, binning or smoothing of the
#' values). The conventions follow the MATLAB workflow of Volponi (2024): samples before
#' the arrival time are set to zero and later deleted, and the fit window ends at `t_end`.
#'
#' * `noise_sigma()` — robust noise SD (1.4826 × MAD) of the background samples, floored
#'   at the quantisation SD of the sensor (`resolution / sqrt(12)`) for flat traces.
#' * `find_hit_time()` — arrival time: the start of the final contiguous run of samples
#'   above `max(3 sigma, frac * peak)` that ends at the peak. Walking back from the peak is
#'   robust to isolated pre-arrival blips; the peak-relative floor rejects small offset
#'   shifts (0.1–0.8 µg/L) seen right after release, which sit above 3σ but far below the
#'   true leading edge.
#' * `find_tail_cutoff()` — the last time before a running-median-smoothed curve stays
#'   below a threshold (LOD = 3σ or LOQ = 10σ) for the rest of the record (at least
#'   `min_run` samples). `NA` when the curve never stays below the threshold, i.e. when
#'   the record ends before the tail reaches the detection limit.
#'
#' @param x background (pre-release) concentrations for `noise_sigma()`.
#' @param resolution sensor quantisation step in concentration units (e.g. ADC step × gain).
#' @param t,C time and background-corrected concentration vectors, in time order.
#' @param sigma noise SD from `noise_sigma()`. @param frac peak fraction floor (0.005).
#' @param thr threshold (e.g. `3 * sigma` or `10 * sigma`).
#' @param min_run minimum number of trailing samples below `thr`. @param k running-median window.
#' @return `noise_sigma()`: list `sigma`, `sigma_mad`, `sigma_sd`, `source`;
#'   `find_hit_time()`: list `hit_time_s`, `threshold`; `find_tail_cutoff()`: a time or `NA`.
#' @examples
#' t <- seq(-300, 2400, by = 30)
#' C <- 300 * dgamma(pmax(t, 0) / 100, 6, 1) + rnorm(length(t), 0, 0.02)
#' s <- noise_sigma(C[t < 0], resolution = 0.015)
#' find_hit_time(t, C, s$sigma)
#' find_tail_cutoff(t, C, 3 * s$sigma)
#' @name tail_rules
NULL

#' @rdname tail_rules
#' @export
noise_sigma <- function(x, resolution = 0) {
  s_mad <- stats::mad(x)
  s_sd <- stats::sd(x)
  floor <- resolution / sqrt(12)
  list(sigma = max(s_mad, floor), sigma_mad = s_mad, sigma_sd = s_sd,
       source = if (s_mad >= floor) "background MAD" else "quantisation floor")
}

#' @rdname tail_rules
#' @export
find_hit_time <- function(t, C, sigma, frac = 0.005) {
  stopifnot(length(t) == length(C), !is.unsorted(t))
  ipk <- which.max(C)
  thr <- max(3 * sigma, frac * C[ipk])
  i <- ipk
  while (i > 1 && C[i - 1] > thr && t[i - 1] > 0) i <- i - 1
  list(hit_time_s = t[i], threshold = thr)
}

#' @rdname tail_rules
#' @export
find_tail_cutoff <- function(t, C, thr, min_run = 5L, k = 5L) {
  stopifnot(length(t) == length(C), !is.unsorted(t))
  ipk <- which.max(C)
  n <- length(C)
  if (ipk >= n) return(NA_real_)
  Cs <- stats::runmed(C, k, endrule = "keep")
  below <- Cs < thr
  for (i in seq(ipk + 1L, n)) {
    if (all(below[i:n]) && (n - i + 1L) >= min_run) return(t[i - 1L])
  }
  NA_real_
}
