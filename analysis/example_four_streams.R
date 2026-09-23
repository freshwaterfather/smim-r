# analysis/example_four_streams.R
# Standalone worked example: fit the SMIM (Volponi et al. 2025 form: V, D, Lambda, beta
# free; t1 = 1.5 s, t2 = 14 000 s fixed) to four rhodamine WT breakthrough curves
# (streams A-D, one release) and draw the chapter figure.
#
# Run from the repository root after installing the package:
#   remotes::install_github("freshwaterfather/smim-r")
#   source("analysis/example_four_streams.R")
#
# Inputs:  data/example_btc.csv     site, t_s (s since release), C_ugL (background-corrected)
#          data/example_config.csv  site, reach_length_m, discharge_L_s, tracer_mass_mg,
#                                   inject_duration_s, hit_time_s, fit_t_end_s, noise_sigma_ugL
# Outputs: figures/example_four_streams_{linear,loglog}.{pdf,png,tiff}   BTCs and fits
#          figures/example_four_streams_params.{pdf,png,tiff}            fitted parameters per stream
#          analysis/example_four_streams_results.csv

library(smimr)
library(ggplot2)

btc <- read.csv("data/example_btc.csv")
cfg <- read.csv("data/example_config.csv")
dir.create("figures", showWarnings = FALSE)

results <- list(); obs <- list(); fitc <- list()
for (s in cfg$site) {
  ci <- cfg[cfg$site == s, ]
  di <- btc[btc$site == s, ]; di <- di[order(di$t_s), ]

  # 1. Prepare exactly as the MATLAB original does (release onward, negatives -> 0,
  #    pre-arrival -> 0, cut at fit_t_end_s, unit-area normalisation, drop t = 0 and C <= 0)
  prep <- prepare_btc(di$t_s, di$C_ugL, hit_time_s = ci$hit_time_s, t_end_s = ci$fit_t_end_s,
                      Q_L_s = ci$discharge_L_s, mass_ug = ci$tracer_mass_mg * 1000)

  # 2. Multi-start bounded least squares: 60 Latin-hypercube starts + the mid-point,
  #    objective of Volponi (1/C-weighted absolute residuals); use objective = "log" for
  #    log-residuals.
  fit <- smim_fit(prep$tobs, prep$cobs, L = ci$reach_length_m, objective = "volponi",
                  bounds = smim_bounds(fix_t1_t2 = TRUE), starts = 60L, seed = 1L, verbose = FALSE)

  results[[s]] <- data.frame(site = s, n_fitted = length(prep$tobs),
                             mass_recovery_pct = 100 * prep$mass_recovery_fraction,
                             v_m_s = fit$params["v"], D_m2_s = fit$params["D"],
                             Lambda_1_s = fit$params["Lambda"], beta = fit$params["beta"],
                             se_v = fit$se_volponi[1], se_D = fit$se_volponi[2],
                             se_Lambda = fit$se_volponi[3], se_beta = fit$se_volponi[4],
                             resnorm = fit$resnorm, r2_lin = fit$diagnostics$r2_lin,
                             r2_log = fit$diagnostics$r2_log)

  # 3. Curves in observation units for the figure: all observations from the release
  #    onward (fitted ones flagged), model on a fine grid from the first sample interval
  area <- pracma::trapz(prep$t_all, prep$C_all)
  tg <- exp(seq(log(min(diff(di$t_s))), log(max(prep$tobs)), length.out = 400))
  fg <- smim_forward(fit$params, tg, ci$reach_length_m)
  fg <- fg / pracma::trapz(prep$tobs, smim_forward(fit$params, prep$tobs, ci$reach_length_m)) * area
  keep <- di$t_s >= 0
  obs[[s]] <- data.frame(site = s, t_s = di$t_s[keep], C = di$C_ugL[keep],
                         in_window = di$t_s[keep] %in% prep$tobs)
  fitc[[s]] <- data.frame(site = s, t_s = tg, C = fg)
}
results <- do.call(rbind, results)
print(results, digits = 4)
write.csv(results, "analysis/example_four_streams_results.csv", row.names = FALSE)

p <- plot_btc_fits(do.call(rbind, obs), do.call(rbind, fitc), x_breaks = c(300, 1000, 2000), x_logticks = TRUE)  # arrivals are > 300 s here
save_btc_figure(p$linear, "figures/example_four_streams_linear", width = 7.2, height = 3.1)
save_btc_figure(p$loglog, "figures/example_four_streams_loglog", width = 7.2, height = 3.1)

# Fitted parameters per stream (layout and labels of Volponi et al. 2025, Fig. 5)
save_btc_figure(plot_smim_params(results), "figures/example_four_streams_params", width = 5.5, height = 4.6)
