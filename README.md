# smimr — Stochastic Mobile-Immobile Model fits to breakthrough curves, in R

An R port of the **SMIMfit** MATLAB library (Schmidt 2020; Volponi 2024) used in
Volponi et al. (2025) to fit the truncated-power-law Stochastic Mobile-Immobile Model
(SMIM) to conservative-tracer breakthrough curves (BTCs). The port reproduces the
original numerics — the de Hoog inverse Laplace transform, the complex incomplete-gamma
waiting-time transform, the 1/C-weighted objective, Latin-hypercube multi-start bounded
least squares and Jacobian-based standard errors — and its equivalence with the MATLAB
original is checked by script (`validation/`).

It accompanies the SMIM section of Chapter 28 (Solute Dynamics in Streams) of *Methods in
Stream Ecology*, 4th edition, and provides the example figure there.

## Install

```r
remotes::install_github("mjliddick/smim-r")
```

## Quick start

```r
library(smimr)
btc <- read.csv("data/example_btc.csv")      # site, t_s, C_ugL (background-corrected)
cfg <- read.csv("data/example_config.csv")   # per-site reach length, discharge, mass, hit time, window
s   <- "A"; ci <- cfg[cfg$site == s, ]; di <- btc[btc$site == s, ]

prep <- prepare_btc(di$t_s, di$C_ugL, hit_time_s = ci$hit_time_s, t_end_s = ci$fit_t_end_s,
                    Q_L_s = ci$discharge_L_s, mass_ug = ci$tracer_mass_mg * 1000)
fit  <- smim_fit(prep$tobs, prep$cobs, L = ci$reach_length_m,
                 objective = "volponi",            # or "log"
                 bounds = smim_bounds(fix_t1_t2 = TRUE),
                 starts = 60L, seed = 1L)
fit$params        # v (m/s), D (m^2/s), Lambda (1/s), beta, log10 t1, log10 t2
fit$se_volponi    # standard errors as computed by the original Lead_SMIM.m
fit$diagnostics   # linear and log R^2, weighted mean absolute error
```

The full worked example (four streams, figure in linear–linear and log–log scale) is
`analysis/example_four_streams.R`.

## What it does

| Step | Function | MATLAB original |
|---|---|---|
| Clean and normalise a BTC (release onward, negatives → 0, pre-arrival → 0, cut at `t_end`, unit area, drop `t = 0` and `C ≤ 0`) | `prepare_btc()` | `Lead_SMIM.m`, `cNorm.m`, `SMIMfit.m` |
| Forward model: travel-time density from the Laplace-domain SMIM solution | `smim_forward()` (`invlap_dehoog()`, `igamma_complex()`) | `TPLmodel.m`, `invLap_deHoog.m`, `LapSolTPL.m`, `memFuncTPL.m`, `LapPsiFuncTPL.m` |
| Objective: `|C_obs − C_fit| / sqrt(N C_obs)` on the area-normalised model (`"volponi"`), or `log C_obs − log C_fit` (`"log"`) | `objective_volponi()`, `objective_log()` | `objFunctionSMIM.m` (+ new `objFunctionSMIM_log.m`) |
| Bounded trust-region-reflective least squares | `trf_lsq()` | `lsqnonlin` (trust-region-reflective) |
| Multi-start from Latin-hypercube points, best local solution refined | `smim_fit()`, `lhs_bounded()` | `runFit.m`, `lhsdesignbnd.m`, `MultiStart` |
| Standard errors | `smim_standard_errors()` | `Lead_SMIM.m` lines 107–140 |
| Figures | `plot_btc_fits()`, `save_btc_figure()` | — |

**Which objective?** `"volponi"` (the default) reproduces the published method: the 1/C
weighting lets the peak dominate, so velocity and peak shape are matched best and the model
may undershoot a long, low tail. `"log"` fits log-residuals, so the many low-concentration
tail points dominate and the tail-derived parameters (Λ, β) follow the observed tail at some
cost near the peak; use it when the goal is to characterise retention from the tail. On the
example data the two give log-scale R² of ~0.98 vs ~0.998 and differ in β by ~30–40 %.

Parameter order everywhere is `c(v, D, Lambda, beta, log10(t1), log10(t2))` (the order the
MATLAB code uses; the comment in `Lead_SMIM.m` lists a different one). By default
`smim_bounds(fix_t1_t2 = TRUE)` fixes `t1 = 1.5 s` and `t2 = 14 000 s` as in
Volponi et al. (2025) and fits the remaining four.

## Validation against MATLAB

`validation/` contains the MATLAB adapter used to run the unmodified original on the
example data (`validation/matlab/`), the exported MATLAB results (`validation/matlab_out/`),
the log of every adaptation (`validation/ADAPTATION_LOG.md`) and the equivalence report
(`validation/VALIDATION_REPORT.md`). The report explains an important property of the
method: the de Hoog inversion at the original's settings is numerically noisy near the
fitted optimum (≈1e-7 of the peak for a 1e-12 relative parameter change), so fitted
parameters are only defined to within a noise-limited basin; the report quantifies this and
states the criteria used.

## Citing

Please cite the original library and paper as well as this port (`CITATION.cff`):

- Volponi, S. N. (2024). SMIMfit library [MATLAB], SMIMfit_2.0. Zenodo. https://doi.org/10.5281/zenodo.11147719 (CC BY 4.0)
- Schmidt, M. J. (2020). SMIMfit library [MATLAB]. https://github.com/mjs271/SMIMfit
- Volponi, S. N., Tank, J. L., Vincent, A. E. S., Snyder, E. D., Pruitt, A. N., & Bolster, D. (2025). Biofilm development, senescence, and benthic substrate influence hyporheic transport in streams. *Journal of Geophysical Research: Biogeosciences*, 130, e2024JG008225. https://doi.org/10.1029/2024JG008225
- Kelly, J. F., et al. (2017). FracFit. *Water Resources Research*, 53, 2559–2567 (the objective's weights).
- de Hoog, F. R., Knight, J. H., & Stokes, A. N. (1982). *SIAM J. Sci. Stat. Comput.*, 3, 357–366 (the inversion).

## License

MIT (see `LICENSE`). The original MATLAB SMIMfit is distributed under CC BY 4.0 on Zenodo;
`lhsdesignbnd.m` (Rik Blok, 2014) is BSD-2. This port is a derivative work with attribution
as required; the MATLAB sources themselves are not redistributed here.
