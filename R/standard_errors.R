#' Parameter standard errors from the final Jacobian
#'
#' Two estimators are provided.
#' * `"volponi"` reproduces `Lead_SMIM.m` lines 107–140 exactly: the error variance is
#'   `sum((cobs - ccfit)^2) / (n - p)` where `ccfit` is the **raw** (not area-normalised)
#'   forward model at the solution, and the covariance is that variance times
#'   `solve(t(J) %*% J)` with `J` the Jacobian of the **weighted** residuals returned by
#'   the optimiser. This mixes an unweighted variance with a weighted Jacobian; it is kept
#'   because it is what the original reports.
#' * `"consistent"` uses the optimiser's own objective: `resnorm / (n - p)` times
#'   `solve(t(J) %*% J)`, the usual Gauss–Newton approximation for the residuals that were
#'   actually minimised.
#' Only the free parameters' Jacobian columns are used (fixed parameters have zero
#' columns, which would make `t(J) %*% J` singular); fixed parameters get `NA`.
#'
#' @param J Jacobian (n x 6) at the solution. @param resid residual vector at the solution.
#' @param cobs,ccfit_raw observed normalised concentrations and the raw forward model.
#' @param free indices of the free parameters.
#' @param method `"volponi"` or `"consistent"`.
#' @return numeric length 6 (NA for fixed parameters).
#' @export
smim_standard_errors <- function(J, resid, cobs, ccfit_raw, free, method = c("volponi", "consistent")) {
  method <- match.arg(method)
  n <- length(cobs); p <- length(free)
  Jf <- J[, free, drop = FALSE]
  cev <- switch(method,
    volponi = sum((cobs - ccfit_raw)^2) / (n - p),
    consistent = sum(resid^2) / (n - p))
  Cov <- cev * solve(crossprod(Jf))
  se <- rep(NA_real_, ncol(J)); se[free] <- sqrt(diag(Cov))
  se
}

#' Goodness-of-fit diagnostics on the fitted points
#'
#' Linear and log-scale coefficients of determination between the observed normalised
#' curve and the area-normalised model curve (the quantity the objective compares), plus
#' the "WMAE" of `Lead_SMIM.m` (the sum of the weighted absolute residuals).
#' @param cobs,ccfit_norm observed and area-normalised fitted concentrations.
#' @param resid residual vector of the objective used.
#' @return list `r2_lin`, `r2_log`, `wmae`.
#' @export
smim_diagnostics <- function(cobs, ccfit_norm, resid) {
  lo <- log(cobs); lf <- log(pmax(ccfit_norm, .Machine$double.xmin))
  list(r2_lin = 1 - sum((cobs - ccfit_norm)^2) / sum((cobs - mean(cobs))^2),
       r2_log = 1 - sum((lo - lf)^2) / sum((lo - mean(lo))^2),
       wmae = sum(resid))
}
