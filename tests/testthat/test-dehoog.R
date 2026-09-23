# de Hoog inversion on transform pairs with closed-form inverses, at the same settings as
# the MATLAB original (tol = 1e-10, M = 25). Accuracy is limited by the method, not the port.

test_that("invlap_dehoog inverts 1/(s+a) -> exp(-a t) across decades", {
  a <- 0.003
  t <- c(seq(30, 990, by = 30), seq(1000, 2800, by = 30))   # two decade groups like our BTCs
  f <- invlap_dehoog(function(s) 1 / (s + a), t)
  expect_equal(f, exp(-a * t), tolerance = 1e-8)
})

test_that("invlap_dehoog inverts 1/s^2 -> t and 1/sqrt(s) -> 1/sqrt(pi t)", {
  t <- seq(100, 2500, by = 50)
  expect_equal(invlap_dehoog(function(s) 1 / s^2, t), t, tolerance = 1e-8)
  expect_equal(invlap_dehoog(function(s) 1 / sqrt(s), t), 1 / sqrt(pi * t), tolerance = 1e-8)
})

test_that("invlap_dehoog returns groups in MATLAB's cell2mat order", {
  t <- c(500, 1500, 50)                    # unsorted: decades 2, 3, 1
  f <- invlap_dehoog(function(s) 1 / s^2, t)
  expect_equal(f, c(50, 500, 1500), tolerance = 1e-8)   # grouped by decade, ascending decade
})

test_that("smim_forward gives a unit-area travel-time density for a benign parameter set", {
  p <- c(0.11, 0.02, 0.03, 1.0, log10(1.5), log10(14000))
  t <- seq(30, 20000, by = 10)
  y <- smim_forward(p, t, L = 60)
  expect_true(all(is.finite(y)))
  expect_gt(pracma::trapz(t, y), 0.9)     # most mass within 20000 s
  expect_lt(pracma::trapz(t, y), 1.02)
})
