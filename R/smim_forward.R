#' Laplace-domain truncated-power-law waiting-time transform
#'
#' Port of `LapPsiFuncTPL.m` (case `'TPL'`). With `tau = T1/T2`, `x = tau + T1 s`,
#' \deqn{\tilde\varphi(s) = [\tau^\beta e^\tau \Gamma(-\beta,\tau)]^{-1}\, x^\beta e^x \Gamma(-\beta, x).}
#' MATLAB compares the complex `x > 1` by real part only; the same rule is used here. The
#' `Re(x) > 1` Gauss–Laguerre branch of the original is reproduced for completeness, but
#' note that the MATLAB original errors in that branch (scalar `denom` indexed by a
#' vector mask), so it can never have contributed to a published fit; a warning is issued.
#'
#' @param s complex vector of Laplace variables.
#' @param param numeric `c(beta, log10(T1), log10(T2))` — the exponentiation is done here,
#'   as in the original.
#' @return complex vector.
#' @keywords internal
lap_psi_tpl <- function(s, param) {
  beta <- param[1]
  T1 <- 10^param[2]
  T2 <- 10^param[3]
  invTau2 <- T1 / T2
  x <- invTau2 + T1 * s
  denom <- (invTau2^beta * exp(invTau2) * igamma_complex(-beta, invTau2))^(-1)
  if (any(Re(x) > 1)) {
    warning("lap_psi_tpl: Re(x) > 1 branch reached; the MATLAB original errors here")
    idx <- Re(x) <= 1
    out <- complex(length(x))
    out[idx] <- denom * x[idx]^beta * exp(x[idx]) * igamma_complex(-beta, x[idx])
    gl <- gauss_laguerre(8L, 0)
    xGL <- x[!idx]
    sum_int <- complex(length(xGL))
    for (i in seq_len(8L)) sum_int <- sum_int + gl$w[i] * (1 + gl$x[i] / xGL)^(-1 - beta)
    out[!idx] <- denom / xGL * sum_int
    return(out)
  }
  denom * x^beta * exp(x) * igamma_complex(-beta, x)
}

#' Laplace-domain exponential (ADE) waiting-time transform, `LapPsiFuncTPL.m` case 'ADE'
#' @keywords internal
lap_psi_ade <- function(s) {
  tchar <- 1
  1 / (1 + tchar * s)
}

#' Gauss–Laguerre nodes and weights (port of `GaussLaguerre.m`, Van Damme 2010)
#' @param n degree; @param alpha generalised parameter (0 = simple Laguerre)
#' @return list with `x` (nodes) and `w` (weights)
#' @keywords internal
gauss_laguerre <- function(n, alpha = 0) {
  i <- seq_len(n)
  a <- (2 * i - 1) + alpha
  b <- sqrt(i[1:(n - 1)] * ((1:(n - 1)) + alpha))
  CM <- diag(a)
  CM[cbind(1:(n - 1), 2:n)] <- b
  CM[cbind(2:n, 1:(n - 1))] <- b
  ev <- eigen(CM, symmetric = TRUE)
  ord <- order(ev$values)
  x <- ev$values[ord]
  V <- ev$vectors[, ord, drop = FALSE]
  w <- gamma(alpha + 1) * V[1, ]^2
  list(x = x, w = as.numeric(w))
}

#' Laplace-domain memory function (port of `memFuncTPL.m`)
#'
#' `options` mirrors the MATLAB cell: `list(c(lambda, tau, W, k_surf, k_sub), "TPL", c(beta, logT1, logT2))`.
#' @keywords internal
mem_func_tpl <- function(s, options) {
  if (length(options) >= 2L) {
    sorp <- options[[1]]
    lambda <- sorp[1]; tau <- sorp[2]; W <- sorp[3]
    if (length(sorp) > 3L) { k_surf <- sorp[4]; k_sub <- sorp[5] } else { k_surf <- 0; k_sub <- 0 }
    lapPsi <- lap_psi_tpl(s + k_sub, options[[3]])
    phi <- W * lapPsi + (1 - W) / (s + k_sub) / tau
    sVec1 <- s + lambda * (1 - phi)
    tchar <- 1
    psi <- lap_psi_ade(sVec1 + k_surf)
    tchar * psi / (1 - psi) * s
  } else {
    1
  }
}

#' Laplace-domain SMIM solution (port of `LapSolTPL.m`, inlet/outlet BC both 'none')
#'
#' @param s complex vector; @param transParams `c(vNorm, DNorm)`;
#' @param options list: `outletBCtype, inletBCfunc, xSample, inletBCtype, sorption params, "TPL", TPL params`.
#' @keywords internal
lap_sol_tpl <- function(s, transParams, options) {
  vNorm <- transParams[1]; DNorm <- transParams[2]
  outletBCtype <- options[[1]]; inletBCfunc <- options[[2]]; xSample <- options[[3]]; inletBCtype <- options[[4]]
  if (vNorm < 0) stop(sprintf("vNorm = %6.5e is negative", vNorm))
  if (DNorm < 0) stop(sprintf("DNorm = %6.5e is negative", DNorm))
  memFunc <- mem_func_tpl(s, options[5:7])
  if (identical(inletBCfunc, "pulse")) {
    LapFuncInletBC <- 1
  } else {
    # user-defined injection: MATLAB eval()s a string in `u`; here a function of s
    LapFuncInletBC <- inletBCfunc(s)
  }
  # expTerm is computed but unused in the original for the 'none' outlet (kept for parity)
  if (!identical(inletBCtype, "none") || !identical(outletBCtype, "none")) stop("Inlet and outlet BCs must both be 'none'")
  const <- vNorm^2 + 4 * s * DNorm / memFunc
  num <- 1 / sqrt(const) * exp((vNorm - sqrt(const)) / 2 / DNorm * xSample)
  denom <- 1
  LapFuncInletBC * num / denom
}

#' SMIM forward model (port of `TPLmodel.m`, conservative case)
#'
#' Returns the modelled breakthrough curve at the observation times, i.e. the travel-time
#' density \eqn{C(t) = v_N\, \mathcal{L}^{-1}\{\tilde C\}(t)} with \eqn{v_N = v/L},
#' \eqn{D_N = D/L^2}. Output is **not** area-normalised; the objective functions do that.
#'
#' @param params numeric length 6: `c(v, D, Lambda, beta, log10(T1), log10(T2))`
#'   (units m/s, m^2/s, 1/s, -, log10 s, log10 s). Note the order; the comment in
#'   `Lead_SMIM.m` line 13 lists it differently from what the code uses.
#' @param tobs observation times (s), all `> 0`.
#' @param L reach length (m).
#' @param inject_duration injection duration (s); `0` = instantaneous pulse.
#' @return numeric vector, same length as `tobs`.
#' @export
smim_forward <- function(params, tobs, L, inject_duration = 0) {
  stopifnot(length(params) == 6L, all(tobs > 0))
  v <- params[1]; D <- params[2]; Lambda <- params[3]; paramsTPL <- params[4:6]
  k_surf <- 0; k_sub <- 0
  vnorm <- v / L
  Dnorm <- D / L^2
  inletBCfunc <- if (inject_duration > 0) function(u) (1 - exp(-inject_duration * u)) / u else "pulse"
  options <- list("none", inletBCfunc, 1, "none", c(Lambda, 1, 1, k_surf, k_sub), "TPL", paramsTPL)
  fun <- function(s) lap_sol_tpl(s, c(vnorm, Dnorm), options)
  c_out <- invlap_dehoog(fun, tobs)
  if (inject_duration > 0) c_out * vnorm / inject_duration else c_out * vnorm
}
