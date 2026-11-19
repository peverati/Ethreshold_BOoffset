# thresholds_Lv: exact critical masses of rotationally excited levels of {X+, X+, e-}

*Last updated: 2026-11-16 08:02*

**STATUS: DONE (run 4 finished 2026-11-15; results in `RESULT.md`, the source of truth).** Set up 2026-11-15 (Paper 1, items 9 and 10); progress and every short test are logged in `PROGRESS.md`.

Results (exact M_c; certified upper bound in RESULT.md):
* J1v1: 17.0942819(10), quote 17.094282(1) (frozen prediction 17.1072 +- 0.0295, exact - pred -0.0129, PASS)
* L3A: 21.376199(1) (21.3818, -0.0056, PASS)
* L4S: 33.914349(2) (33.9175, -0.0032, PASS)
* L1S: 813.54 +- 0.07, **EXPLORATORY** (slow convergence, binding slope 7.9e-9 Eh per unit M)

**J1v1 caveat:** the binding is linear in M - M_c only asymptotically, with visible curvature (E_b/(M-M_c) rises from 3.9e-4 to 6.3e-4 for x = 0.0005 to 0.1). A quadratic-fit root (17.0970) is NOT reliable; the value rests on the converged fixed-basis crossings. Caveats of the blind test: scale optimization ran 4 min before the freeze (no crossing computed); small-basis values were seen after it; the freeze is author-held and was not committed before the run.

## What and why

M_c is the mass M (in electron masses) of the two identical X+ at which a level becomes bound below
E_th = -M/(2(M+1)), the X + (Xe)(1s) threshold. For natural parity, every L and both exchange symmetries
share this threshold, with relative angular momentum lambda = L.

| target | sector (L, exchange symmetry, root) | BO picture | BO M_c | BO+ad M_c |
|---|---|---|---|---|
| (9) `J1v1` | L=1, 'A', root 1 | 1s sigma_g J=1 v=1 (the (1,1) level of muonic molecules) | 18.494475 | 17.523638 |
| (10) `L3A` | L=3, 'A', root 0 | 1s sigma_g J=3 v=0 | 22.599590 | 21.638056 |
| (10) `L4S` | L=4, 'S', root 0 | 1s sigma_g J=4 v=0 | 35.124626 | 34.161564 |
| (10) `L1S` | L=1, 'S', root 0 | 2p sigma_u J=1 v=0 (long-range well) | 815.5587 | n/a |

Symmetry map: nuclear exchange = electronic g/u parity times (-1)^J. So 1s sigma_g with odd J is 'A' and
with even J is 'S'. 2p sigma_u with J=1 is 'S', and with J=0 it is 'A' (L=0).

Literature and consistency checks:
* J1v1: Frolov & Wardlaw (arXiv:1008.3010; EPJD 63, 339 (2011)) find (J=1, v=1) bound up to m_mu ~ 214.72 m_e, i.e.
  M_c ~ 3670.483/214.718 = 17.0944. This is a basis-limited upper estimate; our converged value is 1.6e-4 lower (docs/literature/frolov_2025_2026_notes.md).
* L3A: ttmu (M = 26.585) has a bound J=3 v=0 level; ppmu (8.880) and ddmu (17.752) do not. So
  17.752 < M_c(J=3) < 26.585. The BO+ad value 21.64 is consistent with this.
* L1S: in BO the 2p sigma_u curve is repulsive at short range. It has a shallow polarization well,
  depth 6.08e-5 Eh (unit mass), at R = 12.6 bohr. BO bound levels appear at:
  * J=0: M = 330.62, 1782.69, 4421.2, ...
  * J=1: M = 815.56, 2867.83, ...
  * J=2: M = 1499.5, ...

  So H2+ (M = 1836) has two J=0 levels and one J=1 level in 2p sigma_u in BO. Literature not checked.
  An independent DVR check confirms the J=1 count at M = 815.7. Earlier runs went only up to M = 10, where the
  L=1 'S' sector is unbound, as expected. In this sector there is no bound level below M ~ 800.

## Method

* **BO brackets** (`bo_brackets.py`, done, 276 s; output in `bo_brackets.txt`/`.json`):
  * 1s sigma_g: zero-energy coupling problem of `paper1/bo_offset/coupling.py`, any L, BO and BO+ad.
    Cross-checked with a `bo_h2plus` count_levels bisection (agreement to 2e-10).
  * 2p sigma_u: U_u(R) from `bo_h2plus.electronic` (state '2pu'), spline plus a -9/(4R^4) tail, shot at zero
    energy. BO only, because the adiabatic correction for 2p sigma_u is not implemented. The 2p sigma_u values are
    accurate to about 1e-5 relative (spline of 210 points).
* **Blind prediction** (`predict.py`, FROZEN 2026-11-15 20:17): parameter-free fold with m_v(R), m_r(R).
  This is the same code path as the v = 4 blind test, which passed. Results are in `PREDICTION_J1v1.md` (sha256
  `551b8f2899a56b9056eaac6a42ea8d6b48e78eb22fcfc6235d156ddb5853e887`):
  * J1v1 = **17.1072 +/- 0.0295** (primary). The empirical variant is 17.0991.
  * J3v0 = 21.3818 +/- 0.0295 (secondary).
  * J4v0 = 33.9175 +/- 0.0295 (secondary).

  The error bar is the largest L > 0 calibration residual (L=1 v=0: +0.0295). Do NOT rerun with `--freeze`.
* **Exact** (`run_exact.py`). This imitates paper1/rotational (rotational_critical.py) with the following:
  * Fixed basis, so the matrices are mass-independent. Brent root of E_k(M) - E_th(M). Each crossing M* is a
    rigorous upper bound to M_c, up to rounding.
  * Two compact total-degree scale sets: the root-0 and root-k optima at N = 12 for J1v1. For L3A/L4S:
    root 0, plus the same with 0.6 alpha_X. For L1S: fixed (1,1,1) and (0.5,0.5,1).
  * Nested diffuse sets (atom x X-atom relative motion, plus mirror), as for L = 1, 2. L1S uses larger
    sets (n_d 40-60, alpha_d 0.5 to 0.03), because its level lives at R ~ 13 bohr.
  * Stages:
    1. scales
    2. compact sequence
    3. compact N_D + nested diffuse sets
    4. checks: largest basis at tol 1e-10 and with margin 1e-11 Eh (rounding-safe bound), plus compact
       N_D +/- 2 with all diffuse sets
    5. binding energies at M*(1+x), x = 0.0005 to 0.2 (for L > 0, binding is linear in M - M_c)
* **Analysis** (`analyze.py` -> `RESULT.md`):
  * best bound
  * rounding-safe bound
  * geometric extrapolation in compact N (with all diffuse sets) and diffuse increments
  * threshold-law fit
  * blind-test comparison (sha256 checked), BO/BO+ad offsets

## Solver extension (solvers/threebody)

* `angular.tensor_components(c)` adds the symmetric traceless Cartesian tensor for any L, by harmonic
  projection. It keeps one entry per index multiset, with multinomial weights: 10 entries for L=3, 15 for L=4.
  `overlap/power/kinetic_terms_L` use these weights. For L <= 2 the old flat 3^L lists are kept unchanged,
  so the cached L=1, 2 terms stay valid. The cache is in `threebody/_cache/angular_terms_v1.pkl`, now with
  L=3, 4. Generating the L=4 terms with sympy took about 15 min once.
* New tests in `tests/test_angular.py`:
  * L=3,4 weighted contraction = full 3^L contraction from an independent detracer formula
  * exchange-symmetry union S+A = unsymmetrized, for L=3,4
  * pointwise kinetic integrand check for L=3
  * second-order quadrature check for L=3
* Results: the full suite passes (38 tests, 9.5 min). Validation with `examples/rotational_validation_L34.py`
  (output `validation_h2p_L34.txt`):
  * H2+ (0, J=3): E - ref = +1.2e-11 (extrapolated +1.2e-12) against Tang et al. 2014.
  * H2+ (0, J=4): -3.1e-9 below a J=0..3 Dunham interpolation. This is consistent with the neglected x^4 term,
    and no verified reference is available.

## Cost (measured on short tests, 2026-11-15)

* J1v1, compact N=18 + 3 diffuse, size 4274 (kept 3313), 6 BLAS threads, loaded machine:
  * build 25 s
  * projection 7 s
  * one crossing 39 s (17 eigen-solves)
  * max RSS 1.6 GB (about 11 n^2 doubles)
* L1S, N=10 + 3 diffuse, size 4562 (kept 2471): 42 s for build plus 5 energies.
* The L=3 and L=4 builds are about 3x more expensive per n^2 (2 and 3 channels, more term pairs).
* Estimates from these measurements (n^2 build, n^3 linear algebra), 8 threads, idle machine:

| target | largest basis | est. wall time | peak memory |
|---|---|---|---|
| J1v1 | 7214 (N=20 + 5 diffuse) | ~30 min | ~4.6 GB |
| L3A | 7712 (N=14 + 4 diffuse) | ~45-50 min (+ ~5 min scales) | ~5.2 GB |
| L4S | 7166 (N=12 + 3 diffuse) | ~40 min (+ ~5 min scales) | ~4.5 GB |
| L1S | 8260 (N=12 + 5 diffuse) | ~25-60 min | ~6 GB |
| **total** | | **~2.5 h (allow up to 4-5 h under load)** | **< 6.5 GB** (one job at a time) |

`run_exact.py --target T --dry` prints a more pessimistic per-basis estimate (it extrapolates n^3 from a
small calibration basis). Its TOTAL lines are 3-6 h at 2 threads under load.

## Launch

From the project root:

    nohup caffeinate -dims scripts/paper1_run4.sh > runs_run4.log 2>&1 &

This runs J1v1, then L3A, then L4S, then L1S, then analyze. Each target is resumable: rerun the same
command after a crash. Single target: `THREEBODY_BLAS_THREADS=8 python3 paper1/thresholds_Lv/run_exact.py --target J1v1`.

## Caveats (state honestly in the paper)

* The L3A/L4S compact sizes are smaller than for L = 1, 2 (memory). Their extrapolation errors will be
  larger. The bounds are rigorous regardless.
* L1S: the BO value has no adiabatic correction, and the exact basis choice is new (only a feasibility test).
  Treat the result as exploratory until the compact and diffuse checks converge.
* The J1v1 scale optimization ran 4 min before the prediction was frozen. It computed no crossing (PROGRESS.md).


**Data note (2026-11-16, main session):** `L1S_binding.csv` was computed with the N=10+5 basis, but its x column is measured from M* = 813.6015 (the N=12+5 crossing). Harmless for the exploratory L1S result; keep in mind when plotting.


**Muonic validation (2026-11-16, main session, user request):** `muonic_check.py` -> `muonic_check.csv`, `muonic_check.log`. Direct solves at Frolov 2025's mass ratios reproduce his dd-mu (1,1), dd-mu (2,0) and tt-mu (2,0) binding energies to all printed digits (relative 1-3e-7; compact N = 16, 18, 20 identical to 1e-12). tt-mu (1,1): 8.034498e-3 m.a.u. (45.21 eV).
