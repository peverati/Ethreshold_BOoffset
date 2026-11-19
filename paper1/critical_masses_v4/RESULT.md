# Result: exact critical mass of the v = 4 symmetric S level of {X+, X+, e-}

*Last updated: 2026-11-15 15:19*

Written by `analyze.py` from the output of `run_exact.py`. Fits: `fits_state4.txt`, `fig_state4.png`.

## Blind-test check

* sha256(PREDICTION.md) = 11b446096db3fff83b1b3aa4ee85cf28ff65996a0984fc19e7e5ce2f5bc01e94 matches prediction.sha256: the prediction is unchanged.
* PREDICTION.md frozen 2026-11-14 17:22; first production entry of run_log.txt: 2026-11-15 14:33.
* Short code tests of run_exact.py (2026-11-15, after the freeze) are described in FREEZE_NOTES.md.

## Exact M_c (v = 4)

**M_c = 98.26550 ± 0.00011** (order-3 fit, x <= 0.1; 9 of 12 masses usable).

* fit spread (order 2/3, windows 0.05-0.3, power law): 98.265396 to 98.265515 (± 1.0e-04)
* basis term (largest |dE| of the convergence checks / slope dE_b/dM): 1.5e-06
* threshold-law coefficients sqrt(E_b) = c1 x + c2 x^2 + c3 x^3: [7.67732e-04 8.10474e-05]
* fixed-basis crossings (margin 0, i.e. upper bounds only up to rounding): 98.26522 to 98.26564 for the alpha_s = 1e-4 bases (scatter ~3e-4 = rounding floor 1e-14 Eh / slope); 2 more diffuse members (alpha_s < 1e-4) are at the rounding floor near M_c and give no bound
* rounding-safe upper bound (margin 1e-11 Eh): 98.269623
* fixed-basis slope at M*: 2.9e-11 to 2.4e-10 Eh per mass unit
* The lower side is extrapolated through the threshold law (not rigorous), as for k = 1..3.

Convergence checks (energy change vs the standard basis, and implied |dM|):

| variant | alpha_s | M | dE (Eh) | implied dM |
|---|---|---|---|---|
| tol=1e-10 | 0.001 | 98.2600 | +4.42e-14 | nan |
| tol=1e-10 | 0.001 | 98.3500 | +2.34e-13 | 2.3e-06 |
| tol=1e-10 | 0.001 | 98.7500 | +1.77e-12 | 2.7e-06 |
| tol=1e-10 | 0.0001 | 98.2600 | +3.38e-14 | nan |
| tol=1e-10 | 0.0001 | 98.3500 | +2.31e-12 | 2.3e-05 |
| tol=1e-10 | 0.0001 | 98.7500 | +2.06e-12 | 3.1e-06 |
| compact N=20 | 0.001 | 98.2600 | -1.64e-14 | nan |
| compact N=20 | 0.001 | 98.3500 | +2.35e-14 | 2.3e-07 |
| compact N=20 | 0.001 | 98.7500 | +2.42e-14 | 3.7e-08 |
| compact N=20 | 0.0001 | 98.2600 | +1.60e-14 | nan |
| compact N=20 | 0.0001 | 98.3500 | -1.52e-13 | 1.5e-06 |
| compact N=20 | 0.0001 | 98.7500 | -1.07e-14 | 1.6e-08 |
| compact N=20 all | 0.001 | 98.2600 | +1.48e-14 | nan |
| compact N=20 all | 0.001 | 98.3500 | +1.74e-14 | 1.7e-07 |
| compact N=20 all | 0.001 | 98.7500 | +2.84e-14 | 4.3e-08 |
| compact N=20 all | 0.0001 | 98.2600 | -3.43e-14 | nan |
| compact N=20 all | 0.0001 | 98.3500 | +9.70e-14 | 9.5e-07 |
| compact N=20 all | 0.0001 | 98.7500 | +1.76e-14 | 2.7e-08 |
| diffuse nc=6 | 0.001 | 98.2600 | +1.78e-15 | nan |
| diffuse nc=6 | 0.001 | 98.3500 | -8.22e-15 | 8.0e-08 |
| diffuse nc=6 | 0.001 | 98.7500 | +1.30e-14 | 2.0e-08 |
| diffuse nc=6 | 0.0001 | 98.2600 | +4.04e-14 | nan |
| diffuse nc=6 | 0.0001 | 98.3500 | -6.58e-14 | 6.4e-07 |
| diffuse nc=6 | 0.0001 | 98.7500 | +1.44e-14 | 2.2e-08 |
| diffuse nd=70, nd2=45 | 0.001 | 98.2600 | -2.68e-12 | nan |
| diffuse nd=70, nd2=45 | 0.001 | 98.3500 | +3.91e-14 | 3.8e-07 |
| diffuse nd=70, nd2=45 | 0.001 | 98.7500 | +4.42e-14 | 6.7e-08 |
| diffuse nd=70, nd2=45 | 0.0001 | 98.2600 | -6.58e-14 | nan |
| diffuse nd=70, nd2=45 | 0.0001 | 98.3500 | -4.14e-14 | 4.0e-07 |
| diffuse nd=70, nd2=45 | 0.0001 | 98.7500 | -3.91e-14 | 5.9e-08 |

## Comparison with the blind prediction

| prediction | value | stated ± | exact − predicted | within stated ±? |
|---|---|---|---|---|
| parameter-free, non-perturbative (primary) | 98.2627 | 0.0095 | +0.0028 | **yes** |
| parameter-free, first order | 98.2627 | 0.0119 | +0.0028 | **yes** |
| empirical (a) | 98.2667 | 0.0050 | -0.0012 | **yes** |
| empirical (a') | 98.2708 |  | -0.0053 |  |
| empirical (b) | 98.2652 |  | +0.0003 |  |
| empirical (b') | 98.2653 |  | +0.0002 |  |
| BO+ad (no non-adiabatic correction) | 98.7585 |  | -0.4930 |  |
| BO | 99.7974 |  | -1.5319 |  |

Exact non-adiabatic offset M_c(BO+ad) − M_c(exact) = 0.4930 (predicted 0.4958); residual exact − predicted offset = -0.0028 (v = 1..3: +0.0095, -0.0062, -0.0046).
Offset vs BO: 1.5319 (v = 1..3: 1.556, 1.528, 1.528).

Pass criteria were fixed in PREDICTION.md; 'within' above also allows the exact error bar.
