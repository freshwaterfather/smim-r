#' Fit the SMIM to a prepared breakthrough curve with multi-start bounded least squares
#'
#' R counterpart of `runFit.m`: a bounded trust-region least-squares run ([trf_lsq()])
#' from every start point, selection of the best local solution (lowest sum of squares
#' among runs that converged, else lowest overall — `MultiStart`'s rule), and one final
#' refinement from that point which also yields the Jacobian for standard errors.
#'
#' @param tobs,cobs fitted points from [prepare_btc()]. @param L reach length (m).
#' @param objective `"volponi"` (original) or `"log"`. See [objective_volponi()].
#' @param bounds list from [smim_bounds()].
#' @param starts either a numeric matrix of start points (one per row, 6 columns; e.g.
#'   the matrix exported by the MATLAB driver, for validation), or an integer number of
#'   Latin-hypercube starts to draw with [lhs_bounded()] (the bounds mid-point guess is
#'   appended as the last row, as `runFit.m` does).
#' @param seed seed for [lhs_bounded()] when `starts` is a number.
#' @param inject_duration injection duration (s).
#' @param tol tolerance used for `ftol`, `xtol`, `gtol` (`runFit.m`: `1e-14`).
#' @param max_nfev,max_iter limits per start (`runFit.m`: `10000`).
#' @param verbose print progress.
#' @return list: `params` (best refined), `resnorm`, `residual`, `jacobian`, `ccfit_raw`,
#'   `ccfit_norm`, `se_volponi`, `se_consistent`, `diagnostics`, `per_start`
#'   (data frame: every local solution), `best_start`, `starts`, `bounds`, `objective`.
#' @export
smim_fit <- function(tobs, cobs, L, objective = c("volponi", "log"), bounds = smim_bounds(),
                     starts = 60L, seed = 1L, inject_duration = 0, tol = 1e-14,
                     max_nfev = 10000L, max_iter = 10000L, verbose = TRUE) {
  objective <- match.arg(objective)
  obj <- smim_objective(objective)
  fn <- function(p) obj(p, tobs, cobs, L, inject_duration)

  if (is.numeric(starts) && length(starts) == 1L) {
    X <- lhs_bounded(as.integer(starts), bounds$lower, bounds$upper, seed = seed)
    X <- rbind(X, bounds$guess)
  } else {
    X <- as.matrix(starts); stopifnot(ncol(X) == 6L)
  }
  ns <- nrow(X)
  per <- vector("list", ns)
  for (k in seq_len(ns)) {
    t0 <- proc.time()[["elapsed"]]
    r <- trf_lsq(fn, X[k, ], bounds$lower, bounds$upper, ftol = tol, xtol = tol, gtol = tol,
                 max_nfev = max_nfev, max_iter = max_iter)
    per[[k]] <- data.frame(start = k, t(setNames(r$x, bounds$names)), resnorm = r$resnorm,
                           exitflag = r$exitflag, iterations = r$iterations, nfev = r$nfev,
                           seconds = proc.time()[["elapsed"]] - t0)
    if (verbose) cat(sprintf("start %2d/%d: resnorm %.10g exitflag %d iter %d nfev %d (%.1f s)\n",
                             k, ns, r$resnorm, r$exitflag, r$iterations, r$nfev, per[[k]]$seconds))
  }
  per <- do.call(rbind, per)
  ok <- per$exitflag > 0
  cand <- if (any(ok)) which(ok) else seq_len(ns)
  best <- cand[which.min(per$resnorm[cand])]
  x0 <- as.numeric(per[best, bounds$names])
  fin <- trf_lsq(fn, x0, bounds$lower, bounds$upper, ftol = tol, xtol = tol, gtol = tol,
                 max_nfev = max_nfev, max_iter = max_iter)

  ccfit_raw <- smim_forward(fin$x, tobs, L, inject_duration)
  ccfit_norm <- ccfit_raw / pracma::trapz(tobs, ccfit_raw)
  list(params = setNames(fin$x, bounds$names), resnorm = fin$resnorm, residual = fin$residual,
       jacobian = fin$jacobian, exitflag = fin$exitflag, iterations = fin$iterations, nfev = fin$nfev,
       ccfit_raw = ccfit_raw, ccfit_norm = ccfit_norm,
       se_volponi = smim_standard_errors(fin$jacobian, fin$residual, cobs, ccfit_raw, bounds$free, "volponi"),
       se_consistent = smim_standard_errors(fin$jacobian, fin$residual, cobs, ccfit_raw, bounds$free, "consistent"),
       diagnostics = smim_diagnostics(cobs, ccfit_norm, fin$residual),
       per_start = per, best_start = best, starts = X, bounds = bounds, objective = objective,
       tobs = tobs, cobs = cobs, L = L)
}
