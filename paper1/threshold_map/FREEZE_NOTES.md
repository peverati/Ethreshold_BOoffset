# Freeze notes: threshold map, 14 levels (third pre-specified test)

`PREDICTION_MAP.md` (frozen, SHA-256 recorded, unchanged since) refers to a work log for what was seen before the freeze. This file reproduces those notes.

* **Frozen:** 2026-11-16 09:35 (`predict_map.py --freeze`), SHA-256 of `PREDICTION_MAP.md` 6e494c54232cbadcc3423309babc48fbc8fd7a1a83a5cdb2f3fff8ef99e02fba. Band ±0.0295 = largest |exact − prediction| of the nine levels known before the freeze (level (1,0)); rms 0.0123. The code reproduces the earlier frozen (1,1), (3,0) and (4,0) predictions (asserted).
* **Seen before the freeze:** (i) Frolov's 1992 finite-basis estimates for (1,2), (2,1) and (2,2), and all earlier exact thresholds; (ii) the BO and BO+ad map. No exact three-body calculation of any target had been run, not even a scale optimisation. A non-frozen dry run of `predict_map.py` with identical code preceded the freeze by nine minutes.
* **Observation at the freeze (not used):** all six J > 0 calibration residuals were negative (exact below prediction, −0.003 to −0.030), so a uniform band is conservative for large J and possibly tight for low J and high v.
* **Seen after the freeze:** scale optimisations and short timing tests for (2,1), (1,2) and (5,0), and a pipeline test on a deliberately small basis. These are not results.
