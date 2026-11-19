# Blind prediction: exact critical mass of the v = 4 symmetric S level of {X+, X+, e-}

*Last updated: 2026-11-14 17:22*

**Frozen 2026-11-14 17:22 (local time), before any exact v = 4 calculation was run.** The sha256 of this file is in
`prediction.sha256`; `analyze.py` checks it. This file must not be edited after the exact run starts.
Produced by `predict.py` (full output: `predict_output.txt`, numbers: `prediction.json`).

M_c = mass M (electron masses) of the two identical X+ at which the 5th symmetric S level (v = 4, root k = 4)
becomes bound below E_th = -M/(2(M+1)).

## Inputs (no new fit)

| level | BO M_c | BO+ad M_c (own asymptote) |
|---|---|---|
| v = 4 | **99.7973934** | **98.7585118** |

* Zero-energy coupling problem of `paper_data/bo_offset/coupling.py` (lam = M + 1/2). Cross-check with
  `bo_h2plus` count_levels bisection (ppw 16/32, Richardson): BO 99.7973934,
  BO+ad 98.7585118. The n <= 4 values reproduce `paper_data/critical_masses/bo_critical_masses.csv`
  (see `predict_output.txt`).
* Mass function m_v(R): `paper_data/bo_offset/massfunc/results_fine.csv` (first-principles, converged to 3e-11).

## Adiabatic piece (BO - BO+ad)

| v | exact (both models solved) | first order |
|---|---|---|
| 1 | 1.01828 | 1.01632 |
| 2 | 1.02222 | 1.02117 |
| 3 | 1.03099 | 1.03044 |
| 4 | 1.03888 | 1.03854 |

## Non-adiabatic offset M_c(BO+ad) - M_c(exact): parameter-free fold, and the v = 1..3 test

| v | M_c(BO+ad) | exact offset | fold, 1st order | fold, non-pert. | residual 1st | residual non-pert. |
|---|---|---|---|---|---|---|
| 1 | 10.0631790 | 0.53755 | 0.52569 | 0.52808 | +0.01186 | +0.00947 |
| 2 | 30.1455552 | 0.50617 | 0.51228 | 0.51236 | -0.00611 | -0.00620 |
| 3 | 59.8156817 | 0.49702 | 0.50169 | 0.50167 | -0.00467 | -0.00465 |
| **4** | 98.7585118 | ? | **0.49584** | **0.49582** | | |

Residual = exact offset - predicted offset. <m_v> weighted for v = 4: 0.4907; fraction of the fold from r > 8 bohr: 0.080.

## PREDICTION

| variant | predicted M_c(exact, v = 4) |
|---|---|
| **parameter-free, non-perturbative fold (primary)** | **98.2627 ± 0.0095** |
| parameter-free, first-order fold | 98.2627 ± 0.0119 |
| empirical (a): primary − extrapolated residual (residual line in 1/lam through v = 2, 3 = -0.0040) | 98.2667 ± 0.005 |
| empirical (a'): primary − least-squares residual line over v = 1..3 (-0.0081) | 98.2708 |
| empirical (b): exact BO+ad offsets of v = 1..3 extrapolated in 1/lam (line through v = 2, 3) | 98.2652 |
| empirical (b'): same, quadratic through v = 1..3 | 98.2653 |

**Uncertainties.**
* Parameter-free: the quoted ± is the largest |residual| of v = 1..3 (non-perturbative 0.0095,
  rms 0.0071; first order 0.0119, rms 0.0082). The v = 2, 3 residuals are
  both negative (mean -0.0054), so if that trend holds the exact M_c lies *above* the parameter-free
  value by about 0.005; the parameter-free number does not include it.
* Empirical (a): ± 0.005 is a judgement (the other empirical variants (a'), (b), (b') all lie within it); it is not a statistical error. Empirical variants are fits to the v = 1..3 data and are NOT parameter-free.
* The BO and BO+ad masses themselves are numerically exact to ~1e-8 and do not contribute.

**Pass/fail criterion (set now).** The parameter-free theory passes if |M_c(exact) - 98.2627| <= 0.0095.
The empirical (a) extrapolation passes if |M_c(exact) - 98.2667| <= 0.005.
