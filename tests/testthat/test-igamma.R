# Equivalence building block: complex upper incomplete gamma vs MATLAB igamma
# Reference file produced by validation/matlab/export_igamma_reference.m (Symbolic Math
# Toolbox igamma on 6080 points: 19 beta values x de Hoog s-grids for 6 time scales + generic z).

ref_path <- testthat::test_path("..", "..", "validation", "matlab_out", "igamma_reference.csv")

test_that("igamma_complex matches MATLAB igamma to 1e-12 for |z| <= 2 (model range is |z| < 0.5)", {
  skip_if_not(file.exists(ref_path), "igamma_reference.csv not present")
  ref <- read.csv(ref_path)
  z <- complex(real = ref$re_z, imaginary = ref$im_z)
  G <- complex(real = ref$re_G, imaginary = ref$im_G)
  rel <- rep(NA_real_, nrow(ref))
  for (b in unique(ref$beta)) {
    i <- ref$beta == b
    rel[i] <- Mod(igamma_complex(-b, z[i]) - G[i]) / Mod(G[i])
  }
  expect_lt(max(rel[Mod(z) <= 2]), 1e-12)
  expect_lt(max(rel[Mod(z) <= 20]), 1e-10)   # generic points up to z = 5 (series cancellation)
})

test_that("integer and near-integer orders agree (beta = 1 limit)", {
  z <- c(0.001 + 0.05i, 0.2 + 0.3i, 0.5)
  a1 <- igamma_complex(-1, z)
  a2 <- igamma_complex(-1 + 1e-9, z)
  a3 <- igamma_complex(-1 - 1e-9, z)
  expect_equal(a1, a2, tolerance = 1e-7)
  expect_equal(a1, a3, tolerance = 1e-7)
})

test_that("Gamma(0, x) equals the exponential integral E1(x) for real x", {
  x <- c(0.01, 0.1, 0.5, 1, 2)
  expect_equal(Re(igamma_complex(0, x)), pracma::expint_E1(x), tolerance = 1e-13)
  expect_equal(Im(igamma_complex(0, x)), rep(0, length(x)), tolerance = 1e-15)
})

test_that("recurrence Gamma(a+1,z) = a Gamma(a,z) + z^a e^-z holds", {
  z <- c(0.02 + 0.1i, 0.3 - 0.2i, 0.7)
  for (a in c(-1.2, -1.5, -1.78, -1.9, -1.0001)) {   # a + 1 must stay <= 0
    lhs <- igamma_complex(a + 1, z)
    rhs <- a * igamma_complex(a, z) + z^a * exp(-z)
    expect_equal(lhs, rhs, tolerance = 1e-12)
  }
})
