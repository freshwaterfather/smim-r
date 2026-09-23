# Validation of the R port against the MATLAB original

Original: Volponi (2024) SMIMfit v2.0 (Zenodo 10.5281/zenodo.11147719), run unmodified in
MATLAB R2023b through the adapter in `validation/matlab/` (see `ADAPTATION_LOG.md`).
Port: `smimr` (this repository). All comparisons are scripted; the ones on the example data
run as package tests (`tests/testthat/test-forward-reference.R`, `test-igamma.R`,
`test-dehoog.R`).

## 1. Building blocks

| Component | Reference | Result |
|---|---|---|
| Complex upper incomplete gamma `igamma_complex(-β, z)` | MATLAB `igamma` (Symbolic Math Toolbox) on 6080 points: 19 β values (incl. 0, 1, 2, 1 ± 1e-6) × the de Hoog frequency grids of six time scales, plus generic points | max relative difference 3.7e-15 for \|z\| ≤ 1 (model range \|z\| < 0.5); 8.8e-12 on generic points up to z = 5 |
| de Hoog inversion `invlap_dehoog()` (tol = 1e-10, M = 25, decade grouping) | closed-form pairs 1/(s+a), 1/s², 1/√s | agreement 1e-8 (the method's own accuracy at these settings) |

## 2. Forward model and objectives (equivalence items 1–2)

For each example curve (sites A–D, fitted points as prepared by the original workflow) a
47-row parameter design (40 Latin-hypercube rows over the fit bounds, plus rows at β = 0, 1,
2, 1 ± 1e-6, the bounds mid-point and the bound corners) was evaluated with the MATLAB
`TPLmodel` / `objFunctionSMIM` / `objFunctionSMIM_log` (`validation/matlab/export_forward_reference.m`,
17-significant-digit CSV in `validation/matlab_out/forward_reference/`) and with the R functions.

| site | fitted points | design rows with ≥ 1 % of mass in the window | forward model, max \|ΔC\| / peak | median | objective (original), max rel. diff | objective (log), max rel. diff |
|---|---|---|---|---|---|---|
| A | 71 | 30 | 9.3e-11 | 2.4e-12 | 3.9e-11 | 6.0e-9 |
| B | 55 | 31 | 3.3e-9 | 3.5e-12 | 5.0e-9 | 2.6e-9 |
| C | 73 | 31 | 1.7e-9 | 2.5e-12 | 3.7e-9 | 1.7e-8 |
| D | 76 | 31 | 8.7e-10 | 3.2e-12 | 1.5e-9 | 1.2e-8 |

Tolerances used by the tests: forward model < 1e-6 relative to the curve peak, objectives
< 1e-5 relative, on design rows whose curve carries at least 1 % of the tracer mass inside the
observation window.

**Why these tolerances rather than 1e-8 / 1e-10.** The two languages differ by one unit in the
last place in complex `exp`/`sqrt`. On identical Laplace-domain inputs the two inversions agree
to ≤ 1e-13 on realistic curves, but de Hoog's continued-fraction extrapolation amplifies such
input differences on nearly empty or extremely sharp curves (design rows at the bound corners:
up to 0.5 of the peak, on curves with < 0.1 % of the mass in the window). More importantly,
the inversion at the original's settings is itself noisy near fitted optima: perturbing the
parameters by a relative 1e-12 changes the curve by ~1e-7 of its peak, and the objective by
~1e-6–1e-5 relative. That noise floor, not the port, sets what "equal" can mean.

### Her own example (2795-m reach, `example_Release1`)
The archived parameters re-evaluated in the installed MATLAB return the archived objective to all
digits (toolchain check). On the same data the R and MATLAB `igamma` values agree with a 34-digit
reference to 1e-16 and the Laplace-domain values to 7e-13, and R's inversion fed with MATLAB's
Laplace values reproduces MATLAB's curve to 2e-14 — but the two independently computed curves
differ by 1.2e-4 of the peak, because that sharply peaked curve is ~10⁴ times worse conditioned for
the de Hoog inversion than the example curves here (1e-13 input noise → 2.7e-5 of the peak vs
3.7e-9). Re-fitting her example with the shipped six-parameter settings and 61 seeded starts also
lands in a different local minimum (β ≈ 0.001, log₁₀t₂ ≈ 11.5) than her archived fit (β = 0.28,
log₁₀t₂ = 5.3), which was started from a previous solution; the six-parameter problem is
multimodal, consistent with the paper's choice to fix t₁ and t₂. Details: `ADAPTATION_LOG.md` §1.

## 3. Fitted parameters (equivalence item 3)

Fits were compared on four curves (not all from the published example) with identical start
points, bounds (t₁ = 1.5 s and t₂ = 14 000 s fixed, four free parameters), tolerances (1e-14)
and evaluation limits (10 000): R's objective evaluated at MATLAB's solution agrees with
MATLAB's reported value to 1e-6–2e-5 relative on all four; on three curves every fitted
parameter agrees to within 0.02–0.94 of the standard error reported by the original and the
best objective values agree to 3e-6–1.5e-3; on the fourth both codes settle in a long, shallow
valley of the objective (Λ and β differ by ~20 % at < 2 % difference in objective, and R with
more starts reaches a lower value than MATLAB), which is a property of that curve under the
fixed-t₂ model, not of the port. Because a single start in ~60 reaches the objective floor on
such curves, the "best of 61" is a random draw; the per-start distributions of the two
languages coincide (same floor, same median, same fraction of starts near the floor). The
detailed per-curve tables are kept with the private data.

## 4. Consequence for users

The fitted parameters of this method are defined only to within the numerical noise basin of
the de Hoog inversion (≈ 0.1–1 % of the parameter values on well-behaved curves, more on
curves with a shallow objective valley). Differences between fits smaller than that — in
either language, or between two random seeds — are not meaningful. The standard errors from
the Jacobian do not include this component.
