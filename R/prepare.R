#' Prepare a breakthrough curve for fitting, following Volponi's `Lead_SMIM.m` + `SMIMfit.m`
#'
#' Reproduces, in order, the data handling of the original MATLAB workflow so that the R
#' fit sees exactly the same numbers as the MATLAB fit:
#' 1. keep samples with `t >= 0` (release onward; `Lead_SMIM.m` lines 68–73);
#' 2. set negative concentrations to zero (`Lead_SMIM.m` line 57) — this only affects the
#'    normalisation integral, because those samples are deleted in step 7;
#' 3. set samples before `hit_time` to zero (`Lead_SMIM.m` lines 76–77);
#' 4. drop samples after `t_end` (`Lead_SMIM.m` lines 85–88);
#' 5. normalise to unit area with the trapezoidal rule (`cNorm.m`, which reduces to
#'    `c / trapz(t, c)` whatever `Q` is) and compute the mass-recovery fraction
#'    `Q * trapz(t, c) / mass`;
#' 6. drop the sample at `t == 0` (`SMIMfit.m` lines 127–130);
#' 7. drop every sample with `C <= 0` (`SMIMfit.m` lines 133–135).
#'
#' No other transformation is applied: no smoothing, binning, clamping of the fitted
#' values or imputation. The concentration column must already be gain- and
#' background-corrected.
#'
#' @param t elapsed time since release (s). @param C corrected concentration (µg L⁻¹).
#' @param hit_time_s arrival time (s); samples before it are zeroed (then dropped).
#' @param t_end_s last time retained (s). Default: the last sample.
#' @param Q_L_s discharge (L s⁻¹) and `mass_ug` injected mass (µg), for the recovery
#'   fraction only (they do not affect the normalised curve).
#' @return list with `tobs`, `cobs` (fitted points), `t_all`, `C_all`, `cnorm_all` (the
#'   cleaned window before steps 6–7), `mass_recovery_fraction`, and counts
#'   `n_neg_clamped`, `n_prehit_zeroed`, `n_dropped_le0`.
#' @export
prepare_btc <- function(t, C, hit_time_s = -Inf, t_end_s = max(t), Q_L_s = NA, mass_ug = NA) {
  stopifnot(length(t) == length(C), !is.unsorted(t))
  keep <- t >= 0
  t <- t[keep]; C <- C[keep]
  n_neg <- sum(C < 0)
  C[C < 0] <- 0
  n_prehit <- sum(t < hit_time_s & C > 0)
  C[t < hit_time_s] <- 0
  cut <- t > t_end_s
  t <- t[!cut]; C <- C[!cut]
  area <- pracma::trapz(t, C)
  cnorm <- C / area
  mrf <- if (is.na(Q_L_s) || is.na(mass_ug)) NA_real_ else Q_L_s * area / mass_ug
  tt <- t; cc <- cnorm
  if (tt[1] == 0) { tt <- tt[-1]; cc <- cc[-1] }
  drop <- cc <= 0
  list(tobs = tt[!drop], cobs = cc[!drop], t_all = t, C_all = C, cnorm_all = cnorm,
       mass_recovery_fraction = mrf, n_neg_clamped = n_neg, n_prehit_zeroed = n_prehit,
       n_dropped_le0 = sum(drop))
}

#' Default parameter bounds and start point of Volponi's `Lead_SMIM.m`
#'
#' Order `c(v, D, Lambda, beta, log10(T1), log10(T2))`. With `fix_t1_t2 = TRUE` (the form
#' used in Volponi et al. 2025, §2.4) `log10(T1)` and `log10(T2)` are fixed at
#' `log10(1.5)` and `log10(14000)` by equal lower and upper bounds. The start point is
#' the bounds mid-point, as in `Lead_SMIM.m` line 100.
#' @param fix_t1_t2 logical. @param t1_s,t2_s the fixed tempering times (s).
#' @return list `lower`, `upper`, `guess`, `free` (indices of free parameters), `names`.
#' @export
smim_bounds <- function(fix_t1_t2 = TRUE, t1_s = 1.5, t2_s = 14000) {
  lower <- c(0.05, 2e-10, 1e-4, 0, 0, 0.5)
  upper <- c(1.5, 2.3, 0.5, 2, 0.2, 12)
  if (fix_t1_t2) {
    lower[5] <- upper[5] <- log10(t1_s)
    lower[6] <- upper[6] <- log10(t2_s)
  }
  list(lower = lower, upper = upper, guess = (lower + upper) / 2,
       free = which(lower < upper), names = c("v", "D", "Lambda", "beta", "logT1", "logT2"))
}
