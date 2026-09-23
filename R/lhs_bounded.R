#' Latin-hypercube start points within bounds (port of `lhsdesignbnd.m`)
#'
#' `lhsdesignbnd.m` (Rik Blok, 2014, BSD-2) draws a unit Latin hypercube with MATLAB's
#' `lhsdesign` (default `'criterion','maximin','iterations',5`) and rescales each column
#' to `[lb, ub]`. Here the unit design comes from `lhs::maximinLHS()` when available,
#' otherwise `lhs::randomLHS()`, otherwise a plain stratified draw. Columns with
#' `lb == ub` are constant. R's and MATLAB's random streams differ, so this function is
#' for **independent** runs; for validation against MATLAB use the exported start-point
#' matrix instead (see [smim_fit()] argument `starts`).
#'
#' @param n number of points. @param lb,ub bounds (length p).
#' @param seed integer seed (`set.seed`), recorded in the attribute `"seed"`.
#' @return `n x p` matrix with attribute `seed`.
#' @export
lhs_bounded <- function(n, lb, ub, seed = 1L) {
  p <- length(lb); stopifnot(length(ub) == p, all(ub >= lb))
  set.seed(seed)
  U <- if (requireNamespace("lhs", quietly = TRUE)) lhs::maximinLHS(n, p) else {
    sapply(seq_len(p), function(j) (sample.int(n) - stats::runif(n)) / n)
  }
  X <- sweep(sweep(U, 2, ub - lb, "*"), 2, lb, "+")
  colnames(X) <- names(lb)
  attr(X, "seed") <- seed
  X
}
