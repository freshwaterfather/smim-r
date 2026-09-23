#' Bounded nonlinear least squares by a trust-region reflective method
#'
#' An R implementation of the trust-region reflective (TRF) algorithm of Branch, Coleman
#' & Li (1999) for bound-constrained least squares, following the structure of SciPy's
#' `least_squares(method = "trf")` (the "exact" SVD-based trust-region subproblem, the
#' reflective step selection and the Coleman–Li scaling). It is the R counterpart of
#' MATLAB's `lsqnonlin` default `'trust-region-reflective'` algorithm used in Volponi's
#' `runFit.m`. The two implementations share the Coleman–Li affine scaling and the
#' reflective step but differ in the subproblem solver (MATLAB: 2-D subspace with
#' preconditioned conjugate gradients) and in the exact form of the stopping rules, so
#' iterates are not bit-identical; equivalence is judged on the converged solutions.
#'
#' Jacobians are formed by central finite differences with MATLAB's default step,
#' `eps^(1/3) * max(|x|, 1)`, kept inside the bounds by falling back to a one-sided
#' difference when a central step would cross a bound.
#'
#' @param fun function returning the residual vector for a parameter vector.
#' @param x0 starting point (moved strictly inside the bounds if necessary).
#' @param lb,ub bounds (finite or infinite). Parameters with `lb == ub` are held fixed.
#' @param ftol,xtol,gtol tolerances on the relative cost change, the step, and the scaled
#'   gradient (MATLAB's `FunctionTolerance`, `StepTolerance`, `OptimalityTolerance`).
#' @param max_nfev maximum number of residual evaluations (counting Jacobian columns),
#'   like MATLAB's `MaxFunctionEvaluations`; `max_iter` like `MaxIterations`.
#' @param fd_step finite-difference base step (`eps^(1/3)` = MATLAB's central default).
#' @param verbose print one line per iteration.
#' @return list: `x`, `resnorm` (sum of squares), `residual`, `jacobian` (central FD at
#'   `x`, same rule as MATLAB), `exitflag` (1 gtol, 2 ftol, 3 xtol, 0 limit reached),
#'   `iterations`, `nfev`, `optimality`, `message`.
#' @references Branch, M. A., Coleman, T. F., & Li, Y. (1999). A subspace, interior, and
#'   conjugate gradient method for large-scale bound-constrained minimization problems.
#'   SIAM J. Sci. Comput. 21, 1–23. Coleman, T. F., & Li, Y. (1996). An interior trust
#'   region approach for nonlinear minimization subject to bounds. SIAM J. Optim. 6, 418–445.
#' @export
trf_lsq <- function(fun, x0, lb = -Inf, ub = Inf, ftol = 1e-8, xtol = 1e-8, gtol = 1e-8,
                    max_nfev = Inf, max_iter = 1000L, fd_step = .Machine$double.eps^(1/3),
                    verbose = FALSE) {
  n_all <- length(x0)
  lb <- rep_len(lb, n_all); ub <- rep_len(ub, n_all)
  stopifnot(all(lb <= ub))
  fixed <- lb == ub
  free <- which(!fixed)
  x_all <- as.numeric(x0); x_all[fixed] <- lb[fixed]
  n <- length(free)
  if (n == 0L) stop("trf_lsq: no free parameters")

  nfev <- 0L
  full_fun <- function(xf) { x <- x_all; x[free] <- xf; nfev <<- nfev + 1L; as.numeric(fun(x)) }
  lbf <- lb[free]; ubf <- ub[free]

  jac <- function(xf, f0) {
    J <- matrix(0, length(f0), n)
    for (i in seq_len(n)) {
      h <- fd_step * max(abs(xf[i]), 1)
      xp <- xf; xm <- xf
      xp[i] <- xf[i] + h; xm[i] <- xf[i] - h
      if (xp[i] > ubf[i]) {            # one-sided backward
        J[, i] <- (f0 - full_fun(xm)) / h
      } else if (xm[i] < lbf[i]) {     # one-sided forward
        J[, i] <- (full_fun(xp) - f0) / h
      } else {
        J[, i] <- (full_fun(xp) - full_fun(xm)) / (2 * h)
      }
    }
    J
  }

  x <- .make_strictly_feasible(x0[free], lbf, ubf, rstep = 1e-10)
  f <- full_fun(x)
  J <- jac(x, f)
  m <- length(f)
  g <- as.numeric(crossprod(J, f))
  cost <- 0.5 * sum(f^2)

  sc <- .cl_scaling(x, g, lbf, ubf)
  Delta <- sqrt(sum((x / sqrt(sc$v))^2))
  if (Delta == 0) Delta <- 1
  alpha <- 0
  status <- NA_integer_
  iteration <- 0L
  step_norm <- 0; actual_reduction <- 0; ratio <- 0

  repeat {
    sc <- .cl_scaling(x, g, lbf, ubf); v <- sc$v; dv <- sc$dv
    g_norm <- max(abs(g * v))
    if (g_norm < gtol) { status <- 1L; break }
    if (verbose) cat(sprintf("iter %4d nfev %5d cost %.10e step %.3e opt %.3e\n", iteration, nfev, cost, step_norm, g_norm))
    if (iteration >= max_iter || nfev >= max_nfev) { status <- 0L; break }

    d <- sqrt(v)
    diag_h <- g * dv
    g_h <- d * g
    J_h <- sweep(J, 2, d, "*")
    # augmented system for the exact trust-region subproblem
    J_aug <- rbind(J_h, diag(sqrt(diag_h), n, n))
    f_aug <- c(f, rep(0, n))
    sv <- svd(J_aug)
    uf <- as.numeric(crossprod(sv$u, f_aug))
    s <- sv$d; V <- sv$v
    theta <- max(0.995, 1 - g_norm)

    actual_reduction <- -1
    while (actual_reduction <= 0 && nfev < max_nfev) {
      tr <- .solve_lsq_trust_region(n, m, uf, s, V, Delta, alpha)
      p_h <- tr$p; alpha <- tr$alpha
      p <- d * p_h
      st <- .select_step(x, J_h, diag_h, g_h, p, p_h, d, Delta, lbf, ubf, theta)
      step <- st$step; step_h <- st$step_h; predicted_reduction <- st$predicted_reduction
      x_new <- .make_strictly_feasible(x + step, lbf, ubf, rstep = 0)
      f_new <- full_fun(x_new)
      step_h_norm <- sqrt(sum(step_h^2))
      if (!all(is.finite(f_new))) { Delta <- 0.25 * step_h_norm; continue }
      cost_new <- 0.5 * sum(f_new^2)
      actual_reduction <- cost - cost_new
      ur <- .update_tr_radius(Delta, actual_reduction, predicted_reduction, step_h_norm, step_h_norm > 0.95 * Delta)
      Delta_new <- ur$Delta; ratio <- ur$ratio
      step_norm <- sqrt(sum(step^2))
      status <- .check_termination(actual_reduction, cost, step_norm, sqrt(sum(x^2)), ratio, ftol, xtol)
      if (!is.na(status)) break
      alpha <- alpha * Delta / Delta_new
      Delta <- Delta_new
    }
    if (actual_reduction > 0) {
      x <- x_new; f <- f_new; cost <- cost_new
      J <- jac(x, f); g <- as.numeric(crossprod(J, f))
    } else {
      step_norm <- 0; actual_reduction <- 0
    }
    iteration <- iteration + 1L
    if (!is.na(status)) break
  }
  if (is.na(status)) status <- 0L

  # final Jacobian at the solution with the same FD rule (for standard errors)
  J_final <- jac(x, f)
  Jf <- matrix(0, m, n_all); Jf[, free] <- J_final
  x_out <- x_all; x_out[free] <- x
  list(x = x_out, resnorm = sum(f^2), residual = f, jacobian = Jf, exitflag = status,
       iterations = iteration, nfev = nfev, optimality = max(abs(g * .cl_scaling(x, g, lbf, ubf)$v)),
       message = c("0" = "iteration/evaluation limit reached", "1" = "gtol satisfied",
                   "2" = "ftol satisfied", "3" = "xtol satisfied")[as.character(status)])
}

# ---- helpers (ports of scipy.optimize._lsq.common / trf) -------------------------------

.make_strictly_feasible <- function(x, lb, ub, rstep = 1e-10) {
  x_new <- x
  active <- .find_active_constraints(x, lb, ub, rstep)
  lower_mask <- active == -1; upper_mask <- active == 1
  if (rstep == 0) {
    x_new[lower_mask] <- .nextafter_up(lb[lower_mask])
    x_new[upper_mask] <- .nextafter_down(ub[upper_mask])
  } else {
    x_new[lower_mask] <- lb[lower_mask] + rstep * pmax(1, abs(lb[lower_mask]))
    x_new[upper_mask] <- ub[upper_mask] - rstep * pmax(1, abs(ub[upper_mask]))
  }
  tight <- (x_new < lb) | (x_new > ub)
  x_new[tight] <- 0.5 * (lb[tight] + ub[tight])
  x_new
}

.nextafter_up   <- function(x) x + abs(x) * .Machine$double.eps + .Machine$double.xmin
.nextafter_down <- function(x) x - abs(x) * .Machine$double.eps - .Machine$double.xmin

.find_active_constraints <- function(x, lb, ub, rtol = 1e-10) {
  active <- integer(length(x))
  if (rtol == 0) {
    active[x <= lb] <- -1L; active[x >= ub] <- 1L
    return(active)
  }
  lower_dist <- x - lb; upper_dist <- ub - x
  lower_threshold <- rtol * pmax(1, abs(lb)); upper_threshold <- rtol * pmax(1, abs(ub))
  lower_active <- is.finite(lb) & (lower_dist <= pmin(upper_dist, lower_threshold))
  active[lower_active] <- -1L
  upper_active <- is.finite(ub) & (upper_dist <= pmin(lower_dist, upper_threshold))
  active[upper_active] <- 1L
  active
}

.cl_scaling <- function(x, g, lb, ub) {
  v <- rep(1, length(x)); dv <- rep(0, length(x))
  mask <- (g < 0) & is.finite(ub)
  v[mask] <- ub[mask] - x[mask]; dv[mask] <- -1
  mask <- (g > 0) & is.finite(lb)
  v[mask] <- x[mask] - lb[mask]; dv[mask] <- 1
  list(v = v, dv = dv)
}

.solve_lsq_trust_region <- function(n, m, uf, s, V, Delta, initial_alpha = 0, rtol = 0.01, max_iter = 10L) {
  phi_and_derivative <- function(alpha, suf, s, Delta) {
    denom <- s^2 + alpha
    p_norm <- sqrt(sum((suf / denom)^2))
    phi <- p_norm - Delta
    phi_prime <- -sum(suf^2 / denom^3) / p_norm
    list(phi = phi, phi_prime = phi_prime)
  }
  suf <- s * uf
  # full-rank check as in scipy
  if (m >= n) {
    threshold <- .Machine$double.eps * m * s[1]
    full_rank <- s[length(s)] > threshold
  } else full_rank <- FALSE
  if (full_rank) {
    p <- -as.numeric(V %*% (uf / s))
    if (sqrt(sum(p^2)) <= Delta) return(list(p = p, alpha = 0, n_iter = 0L))
  }
  alpha_upper <- sqrt(sum(suf^2)) / Delta
  if (full_rank) {
    pd <- phi_and_derivative(0, suf, s, Delta)
    alpha_lower <- -pd$phi / pd$phi_prime
  } else alpha_lower <- 0
  alpha <- initial_alpha
  if (initial_alpha == 0 || !full_rank && initial_alpha == 0) alpha <- max(0.001 * alpha_upper, sqrt(alpha_lower * alpha_upper))
  it <- 0L
  repeat {
    if (alpha < alpha_lower || alpha > alpha_upper) alpha <- max(0.001 * alpha_upper, sqrt(alpha_lower * alpha_upper))
    pd <- phi_and_derivative(alpha, suf, s, Delta)
    if (pd$phi < 0) alpha_upper <- alpha
    ratio <- pd$phi / pd$phi_prime
    alpha_lower <- max(alpha_lower, alpha - ratio)
    alpha <- alpha - (pd$phi + Delta) * ratio / Delta
    it <- it + 1L
    if (abs(pd$phi) < rtol * Delta || it >= max_iter) break
  }
  p <- -as.numeric(V %*% (suf / (s^2 + alpha)))
  p <- p * Delta / sqrt(sum(p^2))
  list(p = p, alpha = alpha, n_iter = it)
}

.intersect_trust_region <- function(x, s, Delta) {
  a <- sum(s^2); b <- sum(x * s); c <- sum(x^2) - Delta^2
  d <- sqrt(b^2 - a * c)
  q <- -(b + sign(b) * d); if (b == 0) q <- -d
  t1 <- q / a; t2 <- c / q
  if (t1 < t2) c(t1, t2) else c(t2, t1)
}

.step_size_to_bound <- function(x, s, lb, ub) {
  non_zero <- s != 0
  s_nz <- s[non_zero]
  steps <- rep(Inf, length(x))
  steps[non_zero] <- pmax((lb - x)[non_zero] / s_nz, (ub - x)[non_zero] / s_nz)
  min_step <- min(steps)
  list(step = min_step, hits = as.integer(sign(s)) * as.integer(steps == min_step))
}

.build_quadratic_1d <- function(J, g, s, diag = NULL, s0 = NULL) {
  v <- as.numeric(J %*% s)
  a <- sum(v^2)
  if (!is.null(diag)) a <- a + sum(s * diag * s)
  a <- 0.5 * a
  b <- sum(g * s)
  if (!is.null(s0)) {
    u <- as.numeric(J %*% s0)
    b <- b + sum(u * v)
    c <- 0.5 * sum(u^2) + sum(g * s0)
    if (!is.null(diag)) { b <- b + sum(s0 * diag * s); c <- c + 0.5 * sum(s0 * diag * s0) }
    return(c(a, b, c))
  }
  c(a, b, 0)
}

.minimize_quadratic_1d <- function(a, b, lb, ub, c = 0) {
  t <- c(lb, ub)
  if (a != 0) { extremum <- -0.5 * b / a; if (lb < extremum && extremum < ub) t <- c(t, extremum) }
  y <- t * (a * t + b) + c
  i <- which.min(y)
  c(t[i], y[i])
}

.evaluate_quadratic <- function(J, g, s, diag = NULL) {
  Js <- as.numeric(J %*% s)
  q <- sum(Js^2)
  if (!is.null(diag)) q <- q + sum(s * diag * s)
  l <- sum(s * g)
  0.5 * q + l
}

.select_step <- function(x, J_h, diag_h, g_h, p, p_h, d, Delta, lb, ub, theta) {
  in_bounds <- function(y) all(y >= lb & y <= ub)
  if (in_bounds(x + p)) {
    p_value <- .evaluate_quadratic(J_h, g_h, p_h, diag = diag_h)
    return(list(step = p, step_h = p_h, predicted_reduction = -p_value))
  }
  sb <- .step_size_to_bound(x, p, lb, ub)
  p_stride <- sb$step; hits <- sb$hits
  r_h <- p_h; r_h[hits != 0] <- -r_h[hits != 0]
  r <- d * r_h
  p <- p * p_stride; p_h <- p_h * p_stride
  x_on_bound <- x + p
  to_tr <- .intersect_trust_region(p_h, r_h, Delta)[2]
  to_bound <- .step_size_to_bound(x_on_bound, r, lb, ub)$step
  r_stride <- min(to_bound, to_tr)
  if (r_stride > 0) {
    r_stride_l <- (1 - theta) * p_stride / r_stride
    if (r_stride == to_bound) r_stride_u <- theta * to_bound else r_stride_u <- to_tr
  } else { r_stride_l <- 0; r_stride_u <- -1 }
  if (r_stride_l <= r_stride_u) {
    q <- .build_quadratic_1d(J_h, g_h, r_h, s0 = p_h, diag = diag_h)
    mq <- .minimize_quadratic_1d(q[1], q[2], r_stride_l, r_stride_u, c = q[3])
    r_stride <- mq[1]; r_value <- mq[2]
    r_h <- p_h + r_h * r_stride
    r <- d * r_h
  } else r_value <- Inf
  p <- p * theta; p_h <- p_h * theta
  p_value <- .evaluate_quadratic(J_h, g_h, p_h, diag = diag_h)
  ag_h <- -g_h; ag <- d * ag_h
  to_tr <- Delta / sqrt(sum(ag_h^2))
  to_bound <- .step_size_to_bound(x, ag, lb, ub)$step
  if (to_bound < to_tr) ag_stride <- theta * to_bound else ag_stride <- to_tr
  q <- .build_quadratic_1d(J_h, g_h, ag_h, diag = diag_h)
  mq <- .minimize_quadratic_1d(q[1], q[2], 0, ag_stride)
  ag_stride <- mq[1]; ag_value <- mq[2]
  ag_h <- ag_h * ag_stride; ag <- ag * ag_stride
  if (p_value < r_value && p_value < ag_value) {
    list(step = p, step_h = p_h, predicted_reduction = -p_value)
  } else if (r_value < p_value && r_value < ag_value) {
    list(step = r, step_h = r_h, predicted_reduction = -r_value)
  } else {
    list(step = ag, step_h = ag_h, predicted_reduction = -ag_value)
  }
}

.update_tr_radius <- function(Delta, actual_reduction, predicted_reduction, step_norm, bound_hit) {
  if (predicted_reduction > 0) ratio <- actual_reduction / predicted_reduction
  else if (predicted_reduction == actual_reduction && actual_reduction == 0) ratio <- 1
  else ratio <- 0
  if (ratio < 0.25) Delta <- 0.25 * step_norm
  else if (ratio > 0.75 && bound_hit) Delta <- Delta * 2
  list(Delta = Delta, ratio = ratio)
}

.check_termination <- function(dF, F, dx_norm, x_norm, ratio, ftol, xtol) {
  ftol_satisfied <- dF < ftol * F && ratio > 0.25
  xtol_satisfied <- dx_norm < xtol * (xtol + x_norm)
  if (ftol_satisfied && xtol_satisfied) return(4L)
  if (ftol_satisfied) return(2L)
  if (xtol_satisfied) return(3L)
  NA_integer_
}
