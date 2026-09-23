test_that("noise_sigma uses MAD and falls back to the quantisation floor", {
  set.seed(1)
  x <- rnorm(20, 0.3, 0.05)
  s <- noise_sigma(x, resolution = 0.015)
  expect_equal(s$sigma, mad(x)); expect_equal(s$source, "background MAD")
  s0 <- noise_sigma(rep(0.3, 20), resolution = 0.015)
  expect_equal(s0$sigma, 0.015 / sqrt(12)); expect_equal(s0$source, "quantisation floor")
})

test_that("find_hit_time ignores pre-arrival blips and offset shifts", {
  t <- seq(-300, 2400, by = 30)
  C <- 3000 * dgamma(pmax(t, 0) / 100, 6, 1)     # peak ~480 near t = 500
  C[t > 0 & t < 300] <- 0.4                        # offset shift after release (> 3 sigma, < 0.5 % of peak)
  C[t == 120] <- 5                                  # isolated blip
  h <- find_hit_time(t, C, sigma = 0.05)
  expect_gte(h$hit_time_s, 300)                      # not fooled by the shift or the blip
  expect_lt(h$hit_time_s, t[which.max(C)])
  expect_equal(h$threshold, 0.005 * max(C))
})

test_that("find_tail_cutoff returns NA when the tail never drops below the threshold", {
  t <- seq(30, 2400, by = 30)
  C <- 300 * dgamma(t / 100, 6, 1) + 1               # 1 unit above baseline for ever
  expect_true(is.na(find_tail_cutoff(t, C, thr = 0.5)))
  C2 <- 300 * dgamma(t / 100, 6, 1)                  # decays to ~0
  te <- find_tail_cutoff(t, C2, thr = 0.5)
  expect_false(is.na(te))
  expect_true(all(C2[t > te + 60] < 0.5))
})
