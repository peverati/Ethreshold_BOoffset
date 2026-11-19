# Blind prediction: exact critical masses of {A+ (M), B+ (rM), e-}

*Last updated: 2026-11-16 09:48*

**Frozen 2026-11-16 09:48 (local time), before any exact crossing search in this family.** sha256 of this file in
`prediction_unequal.sha256`; files chmod 444. Numbers: `prediction_unequal.json`; full model table:
`predict_output.txt`. Produced by `predict.py --freeze`. Exact three-body results seen before the freeze:
see PROGRESS.md (only HD+ validation energies at physical masses and an r = 1 energy comparison; no crossing).

M_c = mass M of the LIGHT nucleus (electron masses; heavy nucleus rM) at which the level becomes bound below
the exact threshold E_th = -mu_Be/2 = -rM/(2(rM+1)) (A+ + (B e)(1s)). Level (J, v) = root v of the
natural-parity L = J sector (no exchange symmetry).

## Method (no fitted parameter; THEORY.md)
Two-channel (1s sigma_g, 2p sigma_u) zero-energy problem with exact g/u couplings (a, K) and the exact
threshold built in, plus second-order non-adiabatic 2x2 mass matrices (Sternheimer, massfunc_gu.csv),
solved non-perturbatively ('2ch+m'). Uncertainty = model-variant spread (u-row masses atomic / frozen
below r = 3) + the symmetric-family residual of the same (L, v) scaled by (1+r)/(2r); floor 0.03.

## PREDICTION (M_c of the light nucleus)

| target | BO | BO+ad (own asymptote) | BO+ad (exact threshold) | 2ch (no masses) | **theory 2ch+m** | null NS | null NT | offset vs BO+ad | change vs previous v |
|---|---|---|---|---|---|---|---|---|---|
| r2_J0v1 | 8.35276 | 7.54956 | 14.88825 | 11.95908 | **11.7133 ± 0.0300** | 7.1464 | 14.4851 | -4.1638 | +nan |
| r2_J0v2 | 23.41750 | 22.60985 | 34.19379 | 29.37793 | **29.0886 ± 0.0300** | 22.2302 | 33.8141 | -6.4788 | -2.3150 |
| r2_J0v3 | 45.67667 | 44.86210 | 61.21013 | 53.78249 | **53.4741 ± 0.0300** | 44.4893 | 60.8374 | -8.6120 | -2.1333 |
| r2_J0v4 | 74.88971 | 74.06909 | 95.88744 | 84.70233 | **84.3858 ± 0.0300** | 73.6993 | 95.5177 | -10.3167 | -1.7047 |
| r2_J1v0 | 4.19479 | 3.42644 | 6.51211 | 5.73142 | **5.5819 ± 0.0368** | 3.1547 | 6.2404 | -2.1555 | +nan |
| r2_J2v0 | 9.59646 | 8.83295 | 11.75262 | 11.30586 | **11.1478 ± 0.0300** | 8.6160 | 11.5357 | -2.3148 | +nan |
| r1.5_J0v1 | 9.25121 | 8.38681 | 14.45777 | 12.12240 | **11.8117 ± 0.0300** | 7.9389 | 14.0098 | -3.4249 | +nan |
| r1.5_J0v2 | 25.98981 | 25.12155 | 34.93253 | 30.96037 | **30.6103 ± 0.0300** | 24.6997 | 34.5107 | -5.4887 | -2.0639 |
| r1.5_J0v3 | 50.72223 | 49.84652 | 63.94402 | 57.67545 | **57.3103 ± 0.0300** | 49.4324 | 63.5299 | -7.4637 | -1.9750 |
| r1.5_J0v4 | 83.18116 | 82.29884 | 101.42930 | 91.80081 | **91.4293 ± 0.0300** | 81.8880 | 101.0185 | -9.1305 | -1.6668 |
| r1.5_J1v0 | 4.63125 | 3.80831 | 6.09976 | 5.59766 | **5.4075 ± 0.0320** | 3.5064 | 5.7978 | -1.5992 | +nan |
| r1.5_J2v0 | 10.63310 | 9.81516 | 11.89434 | 11.64154 | **11.4534 ± 0.0300** | 9.5742 | 11.6533 | -1.6382 | +nan |

* Null NS: BO+ad (own asymptote) minus the exact symmetric-family offset of the same (L, v) (lam units).
  Null NT: the same applied to BO+ad measured from the exact threshold. Null NP (persistence): offset(v) =
  offset(v-1) of this family, evaluated after the exact runs; the theory predicts the level-to-level change
  listed in the last column (offset = M_BOad - M_c).
* **Pass/fail (set now).** A target passes if |M_c(exact) - M_pred| <= unc. The test discriminates where
  |M_pred - null| > unc + (exact error bar): for NS and NT this holds for every level (columns above); for NP
  wherever the predicted change exceeds 0.03 + unc.
* **Secondary (lenient) criterion, set now:** |M_c(exact) - M_pred| <= 0.10. Reason: the floor 0.03 is
  calibrated on the symmetric family only; the O(1/lam) terms of the g/u part (2p pi_u, higher Sigma states as
  explicit channels, potential-like second-order terms near the localization radius) are not calibrated. Even
  0.10 is far smaller than the theory-null differences (>= 0.4 for NT, >= 2.4 for NS, level-to-level offset
  changes 1.7-2.3 for NP).
* Offsets are M_BOad - M_c (negative here: the exact threshold lies below the BO+ad asymptote).
* Primary targets: r = 2, J0v1..J0v4, J1v0, J2v0. Secondary: r = 1.5 (exact run only if time permits).
