#' Residual vectors for the SMIM fit
#'
#' `objective_volponi()` is a line-for-line port of `objFunctionSMIM.m` (Volponi/Schmidt,
#' after Kelly et al. 2017 FracFit): the model curve is normalised to unit area over the
#' observation times with the trapezoidal rule, and the residuals are
#' \eqn{r_i = |C_{obs,i} - C_{fit,i}| / \sqrt{N C_{obs,i}}}, so that a least-squares
#' solver minimises \eqn{\sum (C_{obs} - C_{fit})^2 / (N C_{obs})}. Any `NaN` in the model
#' curve returns a vector of `1e16`. `K_mass` is always 1 in the original and is not
#' exposed.
#'
#' `objective_log()` is the log-residual alternative added for this port (validated
#' against `validation/matlab/objFunctionSMIM_log.m`): \eqn{r_i = \log C_{obs,i} - \log C_{fit,i}}
#' on the same normalised curve, with the same `1e16` penalty when the curve has a `NaN`
#' or a non-positive value.
#'
#' **Which objective?** The two weight the curve differently. The original 1/C weighting
#' lets the high-concentration peak dominate, so velocity and the peak shape are matched
#' best and the model may fall below a long, low tail. Log residuals treat every point in
#' relative terms, so the many low-concentration tail points dominate and the tail-derived
#' parameters (\eqn{\Lambda}, \eqn{eta}) follow the observed tail, at some cost near the
#' peak. Use `"volponi"` to reproduce the published method; use `"log"` when the goal is to
#' characterise retention from the tail. On the example data the two give log-scale
#' \eqn{R^2} of about 0.98 and 0.998 and differ in \eqn{eta} by roughly 30–40 %.
#'
#' @param params numeric length 6, see [smim_forward()].
#' @param tobs,cobs observation times and (area-normalised, strictly positive) concentrations.
#' @param L reach length (m). @param inject_duration injection duration (s).
#' @param forward the forward model, by default [smim_forward()] (injectable for tests).
#' @return numeric residual vector of length `length(cobs)`.
#' @references Kelly, J. F., Bolster, D., Meerschaert, M. M., Drummond, J. D., & Packman,
#'   A. I. (2017). FracFit: A robust parameter estimation tool for fractional calculus
#'   models. Water Resources Research, 53, 2559–2567. \doi{10.1002/2016WR019748}
#' @export
objective_volponi <- function(params, tobs, cobs, L, inject_duration = 0, forward = smim_forward) {
  N <- length(cobs)
  c_fit <- forward(params, tobs, L, inject_duration)
  if (any(is.nan(c_fit))) return(rep(1e16, N))
  c_fit <- c_fit / pracma::trapz(tobs, c_fit)
  wgts <- 1 / sqrt(N * cobs)
  wgts[is.infinite(wgts)] <- 0
  wgts * abs(cobs - c_fit)
}

#' @rdname objective_volponi
#' @export
objective_log <- function(params, tobs, cobs, L, inject_duration = 0, forward = smim_forward) {
  N <- length(cobs)
  c_fit <- forward(params, tobs, L, inject_duration)
  if (any(is.nan(c_fit))) return(rep(1e16, N))
  c_fit <- c_fit / pracma::trapz(tobs, c_fit)
  if (any(c_fit <= 0)) return(rep(1e16, N))
  log(cobs) - log(c_fit)
}

#' Select an objective by name
#' @param objective `"volponi"` (default, the original) or `"log"`.
#' @return the residual function.
#' @export
smim_objective <- function(objective = c("volponi", "log")) {
  switch(match.arg(objective), volponi = objective_volponi, log = objective_log)
}
