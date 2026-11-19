# Freeze notes: (1,1), (3,0) and (4,0) thresholds (second pre-specified test)

`PREDICTION_J1v1.md` (frozen, SHA-256 recorded, unchanged since) refers to a work log for what was seen before the freeze. This file reproduces those notes.

* **Frozen:** 2026-11-15 20:17. Primary prediction (1,1) = 17.1072 ± 0.0295 (parameter-free; empirical variant 17.0991); secondary predictions (3,0) = 21.3818 and (4,0) = 33.9175, same band.
* **Seen before the freeze (caveat).** The scale optimisation of `run_exact.py --target J1v1 --stage 1` (N = 12) ran 20:11–20:13, four minutes before the freeze. Its only exact-energy information: at M = 17.5236 (the BO+ad threshold) root 1 is bound by 2.0e-4 Eh in an N = 12 basis. This says only M_c < 17.52, which every earlier level already implied (exact below BO+ad). No crossing was computed before the freeze. The header's "before any exact calculation" should be read with this caveat.
* **Seen after the freeze:** timing and pipeline tests on deliberately small test bases (bounds 17.0943–17.0944, 21.3762, 33.9144). These are not results.
* The (1,1) scale optimisation of stage 1 (`J1v1_scales.json`) was reused by the production run.
