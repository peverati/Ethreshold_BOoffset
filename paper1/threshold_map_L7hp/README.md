# threshold_map_L7hp: (J, v) = (7, 0) with high-precision (double-double) matrix assembly

*Last updated: 2026-11-17 09:48*

Goal: make the L = 7 threshold of {X+, X+, e-} computable. Run 5 (paper1/threshold_map) skipped J7v0 because
validate_L 7 failed (paper1/threshold_map/validate_L7_diagnostics.md). The frozen prediction for J7v0 is in
paper1/threshold_map/PREDICTION_MAP.md (87.5223 +/- 0.0295) and is never touched here (hash checked before and after).

## Diagnosis (step 1, 2026-11-17)
Precision is lost in the **assembly** of S and H (double-precision summation of ~50-1000 angular terms with powers
up to s^16 per matrix element), not in the 1D Laguerre integrals (exact to 4e-16) or in the eigensolver.
Double assembly leaves elementwise errors of ~1e-15 of the largest element, but S has eigenvalues down to ~1e-15
relative, so double S acquires a spurious negative direction (lam_min/lam_max = -2e-14) and tol 1e-14/1e-15 collapse.
Double-double assembly gives lam_min/lam_max = -4e-16 and monotone variational energies (diag_L7_N8_12.txt).

## Method
* `hp_assembly.py`: `HPCoupledBasis` subclasses `threebody.coupled.CoupledBasis` and only overrides `_assemble`.
  Same-scale 1D factors from exact integer matrices (2D)^T J^p (2D) and exact alpha powers; mixed-scale factors by
  Gauss-Laguerre quadrature in mpmath (80 digits); terms accumulated in vectorized double-double; result rounded to
  double. Projection (canonical orthogonalization) and eigensolver unchanged (double). Parallel over unit pairs
  (fork workers, `HP_NPROC`, default 20).
* Private term cache `_cache/angular_terms_map.pkl` (copy of run 5's at commit 8a41b65) via `tb_hp_local.py`:
  nothing here writes to run 5's cache or to solvers/.
* Term coefficients are stored as floats; 98.6 % are exact 13-smooth rationals. The fix works with the floats.

## Result (2026-11-17)
**M_c(J7v0) = 87.521468 +/- 1e-6**; best upper bound 87.5214689 (N = 12 + 2 diffuse); certified bound (margin 1e-11
at tol 1e-14) 87.5214716. Frozen prediction 87.5223 +/- 0.0295: exact - pred = -0.0008, **PASS** (run 5's blind
test becomes 14/14). Details: J7v0_note.md; full section in paper1/threshold_map/RESULT.md (analyze_map.py reads
data/ here). One check basis (N = 8 + 2 diffuse) had a spurious level from single-pass orthogonalization at tol 1e-14
and is excluded (recheck_guard.py, diag_spurious.txt). A two-pass orthogonalization rerun (data_2pass/) reproduces
the accepted crossings to <= 1e-7 but fails validate_hp, so it is a cross-check only. Run 5 (tol 1e-12) is not
affected: two-pass shifts <= 4.3e-10 (check_run5_orthog.txt).

## Files
| file | purpose |
|---|---|
| hp_assembly.py | double-double assembly (HPCoupledBasis) |
| tb_hp_local.py | private term cache hook |
| diag_precision.py, diag_L7_N8_12.txt | step 1: double vs double-double |
| validate_hp.py, validate_hp.txt | step 3: validate_L with the HP path (L = 7, L = 6 control) |
| run_map_hp.py | explicit copy of run 5's run_map.py with the HP changes (listed in its docstring) + guard (c) |
| recheck_guard.py, data/J7v0_guard.csv | guards applied to every finished row |
| diag_spurious.py, diag_spurious.txt | the N = 8 + 2 diffuse spurious level: single vs two-pass, tol 1e-12..1e-15 |
| check_run5_orthog.py, check_run5_orthog.txt | run 5 unaffected (two-pass shifts <= 4.3e-10) |
| data/, run_J7v0.log | J7v0 production (single pass, tol 1e-14) |
| data_2pass/, run_J7v0_2pass.log | cross-check run with two-pass orthogonalization (not used for results) |
| J7v0_note.md | summary appended to run 5's RESULT.md |
| tests/test_hp.py | tests (dd arithmetic, exact 1D factors, mixed scales, agreement with double at low L, parallel = serial) |

Run tests: `cd paper1/threshold_map_L7hp && python -m pytest -q tests` (14 tests, ~4 s).
