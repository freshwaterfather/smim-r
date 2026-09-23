# Equivalence items 1-2 against the MATLAB original (Volponi's TPLmodel / objFunctionSMIM
# and the log objective file), on the four example curves (sites A-D) and a 47-row
# parameter design (LHS over the fit bounds plus edge rows). Reference files produced by
# validation/matlab/export_forward_reference.m.
#
# Tolerances (see validation/VALIDATION_REPORT.md): the de Hoog inversion amplifies
# 1-ulp differences in complex exp/sqrt, so the forward model is compared relative to each
# curve's peak (< 1e-6) and only on design rows whose curve carries >= 1 % of the tracer
# mass inside the observation window; objectives are compared relative (< 1e-5).

ref_dir <- testthat::test_path("..", "..", "validation", "matlab_out", "forward_reference")
sites <- c("A", "B", "C", "D")

for (s in sites) {
  test_that(sprintf("site %s: forward model and objectives match MATLAB", s), {
    skip_if_not(file.exists(file.path(ref_dir, sprintf("site_%s_design.csv", s))), "reference files not present")
    X  <- as.matrix(read.csv(file.path(ref_dir, sprintf("site_%s_design.csv", s))))
    tc <- read.csv(file.path(ref_dir, sprintf("site_%s_tobs.csv", s)))
    Fm <- as.matrix(read.csv(file.path(ref_dir, sprintf("site_%s_forward.csv", s))))
    ob <- read.csv(file.path(ref_dir, sprintf("site_%s_obj.csv", s)))
    Rv <- as.matrix(read.csv(file.path(ref_dir, sprintf("site_%s_resid_volponi.csv", s))))
    Rl <- as.matrix(read.csv(file.path(ref_dir, sprintf("site_%s_resid_log.csv", s))))

    Fr <- sapply(seq_len(nrow(X)), function(k) smim_forward(X[k, ], tc$tobs, L = 60))
    mass <- apply(Fm, 2, function(y) pracma::trapz(tc$tobs, y))
    realistic <- mass >= 0.01 & mass <= 1.5
    rel_peak <- apply(abs(Fr - Fm), 2, max) / apply(abs(Fm), 2, max)
    expect_lt(max(rel_peak[realistic]), 1e-6)

    rv <- sapply(seq_len(nrow(X)), function(k) objective_volponi(X[k, ], tc$tobs, tc$cobs, L = 60))
    rl <- sapply(seq_len(nrow(X)), function(k) objective_log(X[k, ], tc$tobs, tc$cobs, L = 60))
    expect_equal(colSums(rv^2)[realistic], ob$ss_volponi[realistic], tolerance = 1e-5)
    expect_equal(colSums(rl^2)[realistic], ob$ss_log[realistic], tolerance = 1e-5)
    # residual vectors element-wise (relative to the largest residual of each row)
    scale_v <- apply(abs(Rv), 2, max); scale_l <- apply(abs(Rl), 2, max)
    expect_lt(max((abs(rv - Rv) / rep(scale_v, each = nrow(Rv)))[, realistic]), 1e-5)
    expect_lt(max((abs(rl - Rl) / rep(scale_l, each = nrow(Rl)))[, realistic]), 1e-5)
  })
}

test_that("prepare_btc reproduces the fitted points MATLAB used (site A)", {
  skip_if_not(file.exists(file.path(ref_dir, "site_A_tobs.csv")), "reference files not present")
  btc <- read.csv(testthat::test_path("..", "..", "data", "example_btc.csv"))
  cfg <- read.csv(testthat::test_path("..", "..", "data", "example_config.csv"))
  d <- btc[btc$site == "A", ]; ci <- cfg[cfg$site == "A", ]
  prep <- prepare_btc(d$t_s, d$C_ugL, hit_time_s = ci$hit_time_s, t_end_s = ci$fit_t_end_s)
  tc <- read.csv(file.path(ref_dir, "site_A_tobs.csv"))
  expect_equal(prep$tobs, tc$tobs)
  expect_equal(prep$cobs, tc$cobs, tolerance = 1e-6)   # example file rounds C to 5 decimals
})
