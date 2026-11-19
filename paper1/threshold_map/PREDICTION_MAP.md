# Blind predictions: exact critical masses of the (J, v) threshold map of {X+, X+, e-}

*Last updated: 2026-11-16 09:35*

**Frozen 2026-11-16 09:35 (local time), before any exact three-body calculation of these targets was run.**
sha256 of this file: `prediction_map.sha256`; this file, `prediction_map.json`, `predict_map_output.txt` and the
sha256 file are chmod 444. Produced by `predict_map.py` (BO input `bo_map.json`). Must not be edited after the
exact runs start. Honesty notes on what was seen before the freeze: `PROGRESS.md`.

M_c = mass M (electron masses) of the two identical X+ at which the level becomes bound below
E_th = -M/(2(M+1)). Level (J, v) of the 1s sigma_g manifold = root k = v of the natural-parity L = J sector with
exchange symmetry (-1)^J (odd J 'A', even J 'S').

## Method (no fitted parameter)

BO+ad zero-energy coupling problem (lam = M + 1/2); threshold weights folded with the first-principles mass functions
m_v(R) (radial) and m_r(R) (centrifugal) from `paper1/bo_offset/massfunc/results_fine.csv`, solved
non-perturbatively. Same code path as the two earlier blind tests (v = 4: passed; J1v1, J3v0, J4v0: passed).

## Calibration (nine levels with known exact M_c; used ONLY for the band)

| level | M_c(BO+ad) | offset 1st | offset non-pert | f_cen | predicted | exact | exact - pred |
|---|---|---|---|---|---|---|---|
| (0,1) | 10.063179 | 0.5257 | 0.5281 | 0.000 | 9.5351 | 9.525634 | -0.0095 |
| (0,2) | 30.145555 | 0.5123 | 0.5124 | 0.000 | 29.6332 | 29.63939 | +0.0062 |
| (0,3) | 59.815682 | 0.5017 | 0.5017 | 0.000 | 59.3140 | 59.31866 | +0.0046 |
| (0,4) | 98.758512 | 0.4958 | 0.4958 | 0.000 | 98.2627 | 98.2655 | +0.0028 |
| (1,0) | 4.570565 | 0.3310 | 0.3328 | 0.651 | 4.2378 | 4.2082827 | -0.0295 |
| (1,1) | 17.523638 | 0.4168 | 0.4165 | 0.319 | 17.1072 | 17.094282 | -0.0129 |
| (2,0) | 11.778616 | 0.2779 | 0.2778 | 0.773 | 11.5008 | 11.4894217 | -0.0114 |
| (3,0) | 21.638056 | 0.2562 | 0.2562 | 0.830 | 21.3818 | 21.376199 | -0.0056 |
| (4,0) | 34.161564 | 0.2440 | 0.2441 | 0.864 | 33.9175 | 33.914349 | -0.0032 |

Uncertainty band, fixed now and applied to every target: **+/- 0.0295** = largest |exact - pred| of the nine
(level (1, 0)); rms of the nine residuals 0.0123.

## PREDICTIONS

| level | sector | BO | BO+ad | offset non-pert | f_cen | **parameter-free (non-pert.)** | first order | Frolov 1992 II |
|---|---|---|---|---|---|---|---|---|
| (2,1) | L=2 'S' root 1 | 29.019661 | 28.066887 | 0.3601 | 0.486 | **27.7068 +/- 0.0295** | 27.7064 | 27.715 |
| (1,2) | L=1 'A' root 2 | 40.755428 | 39.766531 | 0.4572 | 0.186 | **39.3093 +/- 0.0295** | 39.3091 | 39.375 |
| (3,1) | L=3 'A' root 1 | 42.355695 | 41.408480 | 0.3275 | 0.588 | **41.0810 +/- 0.0295** | 41.0807 |  |
| (5,0) | L=5 'A' root 0 | 50.318074 | 49.353402 | 0.2361 | 0.887 | **49.1173 +/- 0.0295** | 49.1173 |  |
| (2,2) | L=2 'S' root 2 | 53.968672 | 53.002721 | 0.4107 | 0.328 | **52.5920 +/- 0.0295** | 52.5917 | 53.295 |
| (4,1) | L=4 'S' root 1 | 58.421973 | 57.475700 | 0.3060 | 0.657 | **57.1697 +/- 0.0295** | 57.1695 |  |
| (6,0) | L=6 'S' root 0 | 68.181554 | 67.215424 | 0.2305 | 0.903 | **66.9849 +/- 0.0295** | 66.9850 |  |
| (3,2) | L=3 'A' root 2 | 70.347388 | 69.393851 | 0.3768 | 0.433 | **69.0171 +/- 0.0295** | 69.0168 |  |
| (1,3) | L=1 'A' root 3 | 72.664208 | 71.661196 | 0.4716 | 0.125 | **71.1896 +/- 0.0295** | 71.1895 |  |
| (5,1) | L=5 'A' root 1 | 77.189949 | 76.242851 | 0.2908 | 0.706 | **75.9521 +/- 0.0295** | 75.9519 |  |
| (2,3) | L=2 'S' root 3 | 88.206149 | 87.225441 | 0.4379 | 0.235 | **86.7875 +/- 0.0295** | 86.7874 |  |
| (7,0) | L=7 'A' root 0 | 88.715954 | 87.748556 | 0.2263 | 0.915 | **87.5223 +/- 0.0295** | 87.5223 |  |
| (4,2) | L=4 'S' root 2 | 89.656384 | 88.708986 | 0.3519 | 0.511 | **88.3571 +/- 0.0295** | 88.3569 |  |
| (6,1) | L=6 'S' root 1 | 98.647176 | 97.698595 | 0.2793 | 0.743 | **97.4193 +/- 0.0295** | 97.4192 |  |

**Pass/fail (set now).** Each target passes if |M_c(exact) - prediction| <= 0.0295 (+ the exact error bar).
Summary statistic, also set now: the number of passes out of all targets that reach a converged exact value.
