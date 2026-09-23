#' Upper incomplete gamma function for negative order and complex argument
#'
#' Computes the non-regularised upper incomplete gamma function
#' \eqn{\Gamma(a, z) = \int_z^\infty t^{a-1} e^{-t} dt} for real order
#' \eqn{a \le 0} (typically \eqn{a = -\beta} with \eqn{\beta \in [0, 2]}) and complex
#' argument \eqn{z}, matching MATLAB's Symbolic Math Toolbox `igamma(a, z)`, which
#' Volponi's `LapPsiFuncTPL.m` relies on. Base R has no complex incomplete gamma.
#'
#' The ascending series
#' \deqn{\Gamma(a,z) = \Gamma(a) - \sum_{n \ge 0} \frac{(-1)^n z^{a+n}}{n!\,(a+n)}}
#' is used. At and near non-positive integers \eqn{a = -m} both \eqn{\Gamma(a)} and the
#' \eqn{n = m} term diverge and cancel; that pair is evaluated in closed, cancellation-free
#' form with `expm1`, so the function is accurate for every \eqn{\beta} in the fit bounds,
#' including exactly 0, 1 and 2 (the bounds and the mid-point guess). Principal branches
#' are used for \eqn{z^a} and \eqn{\log z}, as in MATLAB.
#'
#' @param a real scalar order, \eqn{-3 < a \le 0} (checked).
#' @param z complex (or real) vector, \eqn{|z| \le 20}. The series is only used in this
#'   range; the SMIM model needs \eqn{|z| \ll 1}.
#' @return complex vector, same length as `z`.
#' @references de Hoog, F. R., Knight, J. H., & Stokes, A. N. (1982). SIAM J. Sci. Stat.
#'   Comput. 3, 357–366 (context); Temme, N. M. (1996) *Special Functions*, §11.2
#'   (near-integer handling).
#' @export
igamma_complex <- function(a, z) {
  stopifnot(length(a) == 1L, is.numeric(a), a <= 0, a > -3)
  z <- as.complex(z)
  if (any(Mod(z) > 20)) stop("igamma_complex: series form is only used for |z| <= 20")
  if (any(z == 0)) stop("igamma_complex: z = 0 is not supported for a <= 0")

  m   <- round(-a)          # nearest non-positive integer -m
  eps <- a + m              # a = -m + eps, |eps| <= 1/2
  logz <- log(z)            # principal branch
  za   <- exp(a * logz)     # z^a, principal branch

  # ---- regular part of the series: all n != m -----------------------------------------
  # term_n = (-1)^n z^(a+n) / (n! (a+n)); ratio term_{n+1}/term_n = -z (a+n) / ((n+1)(a+n+1))
  S <- complex(length(z))
  zn <- rep(1 + 0i, length(z))          # z^n / n!
  n <- 0L
  repeat {
    if (n != m) {
      term <- ((-1)^n) * zn / (a + n)
      S <- S + term
      if (n > 3 && all(Mod(term) <= 1e-18 * pmax(Mod(S), 1e-300))) break
    }
    n <- n + 1L
    zn <- zn * z / n
    if (n > 600L) stop("igamma_complex: series did not converge")
  }
  S <- za * S

  # ---- singular pair: Gamma(a) - (-1)^m z^(a+m) / (m! (a+m)) ---------------------------
  # Gamma(a) = A(eps)/eps with A(eps) = (-1)^m Gamma(1+eps) / prod_{j=1..m} (j - eps)
  # n=m term = B(eps)/eps with B(eps) = (-1)^m z^eps / m!
  # [A - B]/eps = A0 * [ expm1(eps*D)/eps - expm1(eps*log z)/eps ],  A0 = (-1)^m / m!
  # where D = L(eps)/eps, L = lgamma(1+eps) - sum_j log(1 - eps/j)  (real).
  A0 <- ((-1)^m) / factorial(m)
  D  <- .lgamma_ratio(eps, m)
  pair <- A0 * (.expm1_over_eps(eps * D, eps, D) - .expm1_over_eps(eps * logz, eps, logz))

  pair - S
}

# expm1(x)/eps where x = eps * c ; returns c when eps == 0. Works for complex c.
.expm1_over_eps <- function(x, eps, c) {
  if (eps == 0) return(c + 0i)
  .expm1_complex(x) / eps
}

# complex expm1 with a series for small |x| (R's expm1 is real-only)
.expm1_complex <- function(x) {
  x <- as.complex(x)
  out <- exp(x) - 1                     # fine for |x| >= 1 (no cancellation)
  small <- Mod(x) < 1
  if (any(small)) {
    xs <- x[small]
    # Taylor series sum_{k>=1} x^k / k!, 25 terms: truncation < 1/26! ~ 2.5e-27 at |x| = 1
    s <- xs; term <- xs
    for (k in 2:25) { term <- term * xs / k; s <- s + term }
    out[small] <- s
  }
  out
}

# D(eps) = [lgamma(1+eps) - sum_{j=1}^m log(1 - eps/j)] / eps, accurate as eps -> 0
.lgamma_ratio <- function(eps, m) {
  js <- seq_len(m)
  if (abs(eps) >= 0.25) {
    # direct evaluation loses ~1e-16/|eps| relative accuracy; fine for |eps| >= 0.25
    return((lgamma(1 + eps) - sum(log1p(-eps / js))) / eps)
  }
  # series: log Gamma(1+e) = -gamma_E e + sum_{k>=2} (-1)^k zeta(k) e^k / k
  #         -log(1 - e/j)  = sum_{k>=1} e^k / (k j^k)
  # 30 terms: truncation error < 0.25^30 / 31 ~ 3e-20
  K <- 30L
  gammaE <- 0.57721566490153286
  zk <- .zeta_int(K)
  D <- -gammaE
  for (k in 2:K) D <- D + ((-1)^k) * zk[k] * eps^(k - 1) / k
  for (j in js) for (k in 1:K) D <- D + eps^(k - 1) / (k * j^k)
  D
}

# Riemann zeta at integers 2..K (index k), by direct summation with Euler-Maclaurin tail
.zeta_int <- function(K) {
  out <- rep(NA_real_, K)
  N <- 2000
  n <- seq_len(N)
  for (k in 2:K) {
    s <- sum(n^(-k))
    # tail: integral + Euler–Maclaurin corrections
    tail <- N^(1 - k) / (k - 1) - N^(-k) / 2 + k * N^(-k - 1) / 12
    out[k] <- s + tail
  }
  out
}
