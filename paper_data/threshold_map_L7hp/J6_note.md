# J = 6 thresholds, (6,0) and (6,1): double-double rerun

*Last updated: 2026-11-19 17:26*

**Why.** In run 5 (`paper_data/threshold_map/`, double-precision assembly, tol 1e-12) the J = 6 crossing sequences converged slowly (ratio about 0.85) and were not monotone, and adding diffuse sets *raised* the crossing. This is the double-precision truncation diagnosed for L ≥ 6 (see `README.md` of this folder). Run 5 therefore quoted only one-sided bounds: (6,0) ≤ 66.983757757 and (6,1) ≤ 97.416774043.

**What.** Both levels were recomputed on 2026-11-19 with the double-double assembly of this folder (`run_map_hp.py`, canonical orthogonalization at tol 1e-14, stages 2–5, scales reused from run 5). Launcher `run_J6.sh`, data in `data_J6/` (`J6v*_crossings.csv`, `J6v*_binding.csv`, `J6v*_scales.json`, per-target logs).

**Guard-clean crossings** (margin 0, tol 1e-14; `M_star` in `data_J6/J6v*_crossings.csv`):

* (6,0): compact N = 6, 8, 10: 66.98408545, 66.98373151, 66.98371877; N = 10 + 1 diffuse: **66.98371834** (best accepted bound; run 5: 66.983757757).
* (6,1): compact N = 6, 8, 10, 12: 97.44357292, 97.42011510, 97.41744657, 97.41692694; N = 10 + 1 diffuse: **97.41673744** (best accepted bound; run 5: 97.416774043); N = 8 + 2 diffuse: 97.41684723.

**Rows flagged by the conditioning check only.** Several rows have λ_min/λ_max of the overlap between −1.2e-16 and −3.6e-16 after orthogonalization. The crossing-slope guard (no spurious root) passed for all of them. They are left blank (`nan`) in the crossings files:
(6,0) N = 12: 66.98371731; N = 8, 10, 12 + 2 diffuse: 66.98372871, 66.98371861, 66.98371717.
(6,1) N = 10, 12 + 2 diffuse: 97.41673702, 97.41673227.

**Adopted values.** By the author's decision, these flagged rows were accepted (the flag is at the level of double-precision rounding and the slope guard passed). The two-diffuse sequences N = 8, 10, 12 then converge geometrically (ratios 0.14 and 0.043), and the extrapolation gives

* **M_c(6,0) = 66.983717(1)**
* **M_c(6,1) = 97.416732(2)**

(uncertainty in the last digit, in units of 1e-6). The numerically certified upper bounds (margin 1e-11 Eh) of run 5 remain valid: 66.983798633 and 97.416776689.

**Blind test.** Frozen predictions (`paper_data/threshold_map/prediction_map.json`): 66.9849 and 97.4193 ± 0.0295. Misses (exact − prediction, unrounded prediction 66.984942 and 97.419258): **−0.0012** and **−0.0025**, both PASS. These values are used in the paper (Tables, Fig. 2 and 3) and checked by `paper_data/figures_final/summary_values.py`.
