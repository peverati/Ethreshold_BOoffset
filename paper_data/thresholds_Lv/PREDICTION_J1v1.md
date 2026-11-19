# Blind prediction: exact critical mass of the (J = 1, v = 1) level of {X+, X+, e-}

*Last updated: 2026-11-15 20:17*

**Frozen 2026-11-15 20:17 (local time), before any exact J = 1, v = 1 (or L = 3, 4) three-body calculation was run.**
sha256 of this file: `prediction_J1v1.sha256`; files chmod 444. Produced by `predict.py` (output
`predict_output.txt`, numbers `prediction_J1v1.json`). Must not be edited after the exact runs start.

M_c = mass M (electron masses) of the two identical X+ at which the level becomes bound below
E_th = -M/(2(M+1)). (J=1, v=1) = second root of the L = 1 exchange-antisymmetric natural-parity sector.

## Method (no fitted parameter)
BO+ad zero-energy coupling problem (lam = M + 1/2) for L = 1, second root; threshold weights folded with
m_v(R) (radial) and m_r(R) (centrifugal) from `paper1/bo_offset/massfunc/results_fine.csv`; first order
and non-perturbative. Same code path as the v = 4 blind test (passed: predicted 98.2627 +/- 0.0095, exact
98.26550 +/- 0.00011).

| state | M_c(BO+ad) | offset 1st | offset non-pert | f_cen | <m_r>_w | exact M_c | resid (exact - pred offset) |
|---|---|---|---|---|---|---|---|
| L=0 v=1 | 10.063179 | 0.52569 | 0.52808 | 0.000 | nan | 9.525634 | +0.0095 |
| L=1 v=0 | 4.570565 | 0.33100 | 0.33279 | 0.651 | 0.2215 | 4.2082827 | +0.0295 |
| L=2 v=0 | 11.778616 | 0.27787 | 0.27783 | 0.773 | 0.2107 | 11.4894217 | +0.0114 |
| L=1 v=1 | 17.523638 | 0.41676 | 0.41647 | 0.319 | 0.2401 | ? | |
| L=3 v=0 | 21.638056 | 0.25625 | 0.25622 | 0.830 | 0.2061 | ? | |
| L=4 v=0 | 34.161564 | 0.24404 | 0.24406 | 0.864 | 0.2036 | ? | |

## PREDICTION

| level | BO | BO+ad | **parameter-free (non-pert.)** | first order | empirical |
|---|---|---|---|---|---|
| J1v1 | 18.494475 | 17.523638 | **17.1072 +/- 0.0295** | 17.1069 | 17.0991 |
| J3v0 | 22.599590 | 21.638056 | **21.3818 +/- 0.0295** | 21.3818 | 21.3753 |
| J4v0 | 35.124626 | 34.161564 | **33.9175 +/- 0.0295** | 33.9175 | 33.9133 |

* **Primary: M_c(J=1, v=1) = 17.1072 +/- 0.0295** (parameter-free). J = 3 and J = 4 rows are
  secondary blind predictions made by the same code at the same time.
* Uncertainty = largest |residual| of the two L > 0 calibration levels (L = 1 v = 0: +0.0295;
  L = 2 v = 0: +0.0114). Both residuals are positive (exact offset larger than predicted), i.e. the
  exact M_c tends to lie BELOW the parameter-free value; 'empirical' subtracts a residual
  0.1446/lam (mean of lam*resid of the two calibration levels; a fit, NOT parameter-free).
* Literature cross-check (not used): Frolov & Wardlaw (arXiv:1008.3010), (J=1, v=1) of dd-mu-like system bound up
  to m_mu ~ 214.72 m_e, i.e. M_c ~ 17.09 (basis-limited upper estimate; docs/literature/referee_check_2026-11-15.md).
* Consistency (not used): ttmu (M = 26.585) has a bound J = 3 level, ddmu (17.752) does not: 17.752 < M_c(J=3) < 26.585.

**Pass/fail (set now).** Parameter-free passes if |M_c(exact) - 17.1072| <= 0.0295.
