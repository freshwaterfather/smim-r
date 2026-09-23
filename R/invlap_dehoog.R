#' Numerical inverse Laplace transform (de Hoog, Knight & Stokes 1982)
#'
#' Line-for-line port of Volponi/Schmidt's `invLap_deHoog.m`: the time vector is split
#' into decades (`floor(log10(t))`) and each group is inverted separately with
#' `T = 2 * max(t)`, `gamma = -log(tol) / (2 T)`, `2M + 1` terms of the continued
#' fraction, and the same quotient–difference recurrences and remainder term.
#'
#' @param F function of a complex vector `s` returning the Laplace-domain values.
#' @param t numeric vector of times, all `> 0`.
#' @param tol accuracy parameter (`getTolerance.m`: `1e-10`).
#' @param M number of terms (`getTolerance.m`: `25`).
#' @return numeric vector `f(t)` in the order of `t` **sorted by decade group** — exactly
#'   as MATLAB's `cell2mat` returns it. For a monotone increasing `t` this is the input
#'   order; for other orders the result is permuted the same way MATLAB permutes it.
#' @references de Hoog, F. R., Knight, J. H., & Stokes, A. N. (1982). An improved method
#'   for numerical inversion of Laplace transforms. SIAM J. Sci. Stat. Comput. 3, 357–366.
#'   \doi{10.1137/0903022}
#' @export
invlap_dehoog <- function(F, t, tol = 1e-10, M = 25L) {
  t <- as.numeric(t)
  alpha <- 0
  ns <- 2L * M + 1L

  mags <- floor(log10(t))
  minMag <- min(mags)
  numMags <- max(mags) - minMag + 1L
  out <- vector("list", numMags)

  for (i in seq_len(numMags)) {
    tv <- t[mags == minMag + i - 1L]
    if (length(tv) == 0L) { out[[i]] <- numeric(0); next }   # MATLAB: empty cell, dropped by cell2mat
    nt <- length(tv)
    Tt <- 2 * max(tv)
    gam <- alpha - log(tol) / (2 * Tt)
    sVec <- gam + (1i * pi * (0:(2L * M))) / Tt
    Fs <- as.complex(F(sVec))

    e <- matrix(0 + 0i, ns, M + 1L)
    # MATLAB pre-allocates q as (ns-1) x M but then assigns column M+1 (line 95), which
    # silently grows the matrix; allocate M+1 columns here.
    q <- matrix(0 + 0i, ns - 1L, M + 1L)
    d <- complex(ns)
    A <- matrix(0 + 0i, nt, ns + 1L)
    B <- matrix(1 + 0i, nt, ns + 1L)

    Fs[1L] <- Fs[1L] / 2
    q[, 2L] <- Fs[2:ns] / Fs[1:(ns - 1L)]
    for (r in 2:(M + 1L)) {
      mr <- 2L * (M - r + 1L) + 1L
      e[1:mr, r] <- q[2:(mr + 1L), r] - q[1:mr, r] + e[2:(mr + 1L), r - 1L]
      if (r < M + 1L) {
        rq <- r + 1L
        mr <- 2L * (M - rq + 1L) + 3L
        q[1:(mr - 1L), rq] <- q[2:mr, rq - 1L] * e[2:mr, rq - 1L] / e[1:(mr - 1L), rq - 1L]
      }
    }

    d[1L] <- Fs[1L]
    d[seq(2L, ns - 1L, by = 2L)] <- -q[1L, 2:(M + 1L)]
    d[seq(3L, ns, by = 2L)]      <- -e[1L, 2:(M + 1L)]

    A[, 2L] <- d[1L]
    z <- exp(1i * pi * tv / Tt)
    for (n in 3:ns) {
      A[, n] <- A[, n - 1L] + d[n - 1L] * z * A[, n - 2L]
      B[, n] <- B[, n - 1L] + d[n - 1L] * z * B[, n - 2L]
    }

    h2M <- (1 + (d[ns - 1L] - d[ns]) * z) / 2
    R2M <- h2M * ((1 + d[ns] * z / h2M^2)^0.5 - 1)
    A[, ns + 1L] <- A[, ns] + R2M * A[, ns - 1L]
    B[, ns + 1L] <- B[, ns] + R2M * B[, ns - 1L]

    out[[i]] <- (exp(gam * tv) / Tt) * Re(A[, ns + 1L] / B[, ns + 1L])
  }
  unlist(out)
}
