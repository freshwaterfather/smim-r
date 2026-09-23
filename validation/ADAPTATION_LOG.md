# ADAPTATION_LOG.md — how Volponi's SMIMfit was run on this dataset

Original code: Volponi, S. N. (2024). SMIMfit library [MATLAB], v2.0. Zenodo
<https://doi.org/10.5281/zenodo.11147719> (CC BY 4.0), building on Schmidt (2020)
<https://github.com/mjs271/SMIMfit>. Repository snapshot used: commit `058f87d`.

**None of her files were edited.** Every adaptation lives in new files under
`validation/matlab/`, which put her folder on the MATLAB path read-only. Where her
workflow had to be reproduced rather than called (her driver `Lead_SMIM.m` is a script tied
to an Excel template), the reproduced lines are cited by line number.

Environment: MATLAB R2023b (23.2) with Optimization, Global Optimization, Statistics &
Machine Learning, and Symbolic Math (for `igamma`) Toolboxes. No Parallel Computing
Toolbox; independent `matlab -batch` processes are used instead.

## 1. Reproducing her own example (`repro_her_example.m`)

Purpose: prove the toolchain on her `example_Release1` sheet before touching our data.

| # | Change relative to `Lead_SMIM.m` | Reason |
|---|---|---|
| 1 | Removed line 7, `cd(fileparts(matlab.desktop.editor.getActiveFilename))` | Errors in `matlab -batch` (no editor). Her folder is added with `addpath` instead. |
| 2 | Process only the `example_Release1` sheet (line 27 lists all sheets) | The `Instructions` and `Template` sheets contain no data and would error inside the loop. |
| 3 | `rng(seed,'twister')` immediately before `SMIMfit`; the identical LHS start matrix is regenerated with the same seed and exported to `starts.csv` | Her code never seeds the RNG, so results are not reproducible run to run. `lhsdesign` (inside `lhsdesignbnd`, `runFit.m:51`) is the only RNG consumer, so regenerating after `rng(seed)` yields exactly the points `runFit` draws. |
| 4 | All outputs (`.mat`, results table, curves, PNG) written to `validation/matlab_out/her_example/` | Her script writes into its own folder, which is read-only here. |
| 5 | `saveas(...,'.fig')` replaced by a PNG of the same two panels | No display in batch mode. |
| 6 | `Day` derived from the text dates in column A instead of `xlsread(...,'C:C')` (line 46) | In her own example sheet column C is empty and column A holds text dates; the current `Lead_SMIM.m` therefore errors on its own example (`Day = Day(ConcSize)`: index exceeds 0). Her archived `example_Release1.mat` contains `Day` of the right length and variables (`flagGuess`, `flagMuG`) that do not exist in the current script, i.e. it was produced by an older version of the driver. |

Everything else (bounds, midpoint guess, `cNorm` call, `SMIMfit` call, covariance and
standard-error lines 107–140, results table) is copied verbatim.

Outcome: the run was started but cancelled after 14 of 61 starts when the project scope was
reduced (each start takes ~4 min with six free parameters because `igamma` is evaluated
through the Symbolic engine). The wrapper is provided for anyone who wants to complete it;
`validation/compare_her_example.R` compares its output with her `SMIMResults.xls`.

## 2. Running her code on our BTCs (`load_btc.m`, `run_phase_a.m`)

Her data path is an Excel template read by `Lead_SMIM.m` lines 40–52. Ours is a CSV
produced by a single preprocessing script that applies the sonde gain and the background
correction once, so that MATLAB and R read identical corrected numbers.

| # | Her assumption (file:line) | Our situation | Adaptation |
|---|---|---|---|
| A1 | Excel template, one sheet per site (`Lead_SMIM.m:40-52`) | corrected CSV + fit-window CSV | `load_btc.m` reads them; no corrections applied in MATLAB |
| A2 | Concentration already background-corrected, mg/L (Instructions sheet; `cNorm.m` header) | µg/L (ppb) after `V × gain − pre-release baseline` (done in preprocessing) | units only enter the mass-recovery fraction; the injected mass is passed in µg so `Q [L/s] × C [µg/L] × t [s]` is consistent |
| A3 | Clock times + day numbers, uniform grid rebuilt from `TimeStep` (`Lead_SMIM.m:68-88`) | elapsed seconds since release on a uniform 30-s grid with a sample at t = 0 | elapsed time passed directly; it equals her rebuilt grid |
| A4 | Negatives clamped to 0 before normalisation (`Lead_SMIM.m:57`) | same | reproduced (`run_phase_a.m` step 2). Clamped points are deleted by `SMIMfit.m:133-135` before fitting, so the clamp only affects the normalisation integral |
| A5 | Pre-arrival samples zeroed using a hand-entered `HitTime` (`Lead_SMIM.m:76-77`) | no hand-entered arrival time | `hit_time_s` computed by one rule for all BTCs in preprocessing (start of the final run above max(3σ, 0.5 % of peak) that ends at the peak); reproduced as step 3 |
| A6 | Samples after a hand-entered `CutTime` dropped (`Lead_SMIM.m:85-88`) | manual `fit_t_end_s` for one BTC; otherwise a detection-limit rule or the record end | step 4 uses `t_end` selected by `t_end_mode`; then `SMIMfit(..., 't_end', t(end))` exactly as she does |
| A7 | `SensorLoc` from the sheet (`Lead_SMIM.m:40`, her example 2795 m) | 60 m for all releases | passed as `'L', 60` |
| A8 | `Q` (documentation says m³/s, `cNorm.m` says L/s) and `TracerMass` in mg | Q in L/s, mass 98.9 mg | passed as L/s and 98 900 µg |
| A9 | Six free parameters, bounds `[0.05 2e-10 1e-4 0 0 0.5]`–`[1.5 2.3 0.5 2 0.2 12]`, guess = midpoint (`Lead_SMIM.m:14,17,100`) | the paper (Volponi et al. 2025 §2.4) fixes t1 = 1.5 s, t2 = 14 000 s and fits four parameters | lower = upper = log10(1.5) and log10(14 000) for columns 5–6 (tested: `lsqnonlin` trust-region-reflective accepts equal bounds); other bounds and the midpoint guess unchanged |
| A10 | Standard errors from all six Jacobian columns (`Lead_SMIM.m:110-123`) | two columns are identically zero when t1, t2 are fixed, which makes `JᵀJ` singular | the two fixed columns are dropped and `p = 4`; her formula is otherwise unchanged (`se_volponi.csv`). An internally consistent variant using `resnorm` is exported alongside as a diagnostic (`se_consistent.csv`), not as a replacement |
| A11 | Objective hard-wired to `objFunctionSMIM` via `createModel.m` | a log-residual objective is also required | new file `objFunctionSMIM_log.m` with the same interface; for that objective the driver reproduces `SMIMfit.m:119-142` (t_end cut, drop t = 0, drop C ≤ 0) and calls her `runFit.m` directly |
| A12 | No RNG seed (`runFit.m:51`) | reproducibility required | `rng(seed,'twister')` before the fit; start matrix regenerated and exported (`starts.csv`) |
| A13 | Results kept in the workspace / `.xls` | R must read them | all arrays written as CSV with 17 significant digits |

Not adapted (kept exactly): `restart = 60`, `tol = 1e-14`, `maxiter = 10000`
(`runFit.m:37-39`), de Hoog `tol = 1e-10`, `M = 25` (`getTolerance.m`), the
`objFunctionSMIM` weights and NaN guard, `cNorm` normalisation.

## 3. Known defects in the original that affect the port (not fixed in MATLAB)
- `LapPsiFuncTPL.m:41-63` (`Re(x) > 1` branch): `denom` is a scalar but is indexed with a
  51-element logical mask, so the branch errors whenever it is entered. It is never entered
  for our data (it needs a time group with max t ≲ 9 s). The R port reproduces the branch
  and warns.
- Parameter order in comments (`Lead_SMIM.m:13,16`, documentation PDF) is
  `[V, D, Λ, log10T1, log10T2, β]`; the code uses `[v, D, Λ, β, log10T1, log10T2]`.
- `Lead_SMIM.m:113` computes the error variance from unweighted residuals between `cobs`
  and a `ccfit` that is not area-normalised, while the Jacobian is that of the weighted
  residuals. Ported as-is (`se_volponi`), see A10.
