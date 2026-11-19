# threshold_map: complete (J, v) threshold map of {X+, X+, e-} for M_c <= ~100

*Last updated: 2026-11-16 12:49*

**STATUS: DONE.** Results are in `RESULT.md` (written by `analyze_map.py` after every target); what was seen before the
predictions were frozen is in `FREEZE_NOTES.md`. Level (7,0) and the final (6,0), (6,1) values come from `../threshold_map_L7hp/`.

## What and why

Paper 1 has exact critical masses M_c (mass of the two identical X+, in electron masses, at which a level becomes
bound below E_th = -M/(2(M+1))) for nine levels of the 1s sigma_g manifold: (0,1)-(0,4), (1,0), (1,1), (2,0), (3,0),
(4,0). This folder completes the map: every (J, v) whose BO+ad threshold is <= 101. That gives (i) the full order of
appearance of the levels with converged numbers, (ii) the number of bound levels n(M) exact vs BO vs BO+ad, (iii) a
14-level blind test of the parameter-free mass-function prediction, and (iv) converged values for Frolov's three
low-precision 1992 thresholds (1,2) 39.375, (2,1) 27.715, (2,2) 53.295.

Level (J, v) = root k = v of the natural-parity L = J sector with exchange symmetry (-1)^J (odd J 'A', even J 'S';
nuclear exchange = g/u parity times (-1)^J). The 2p sigma_u levels live in the other symmetry and start at M = 330.6,
so each sector below M ~ 100 contains only the 1s sigma_g ladder: root k is v (checked at every crossing, below).

## Targets (from `bo_map.json`; run order = Frolov's three first, then by BO+ad)

| run | (J,v) | sector | BO | BO+ad | frozen prediction | Frolov 1992 II |
|---|---|---|---|---|---|---|
| 1 | (1,2) | L=1 'A' root 2 | 40.755428 | 39.766531 | 39.3093 +/- 0.0295 | 39.375 (bound at 39.45) |
| 2 | (2,1) | L=2 'S' root 1 | 29.019661 | 28.066887 | 27.7068 +/- 0.0295 | 27.715 (bound at 27.75) |
| 3 | (2,2) | L=2 'S' root 2 | 53.968672 | 53.002721 | 52.5920 +/- 0.0295 | 53.295 (bound at 53.30) |
| 4 | (3,1) | L=3 'A' root 1 | 42.355695 | 41.408480 | 41.0810 +/- 0.0295 | |
| 5 | (5,0) | L=5 'A' root 0 | 50.318074 | 49.353402 | 49.1173 +/- 0.0295 | |
| 6 | (4,1) | L=4 'S' root 1 | 58.421973 | 57.475700 | 57.1697 +/- 0.0295 | |
| 7 | (6,0) | L=6 'S' root 0 | 68.181554 | 67.215424 | 66.9849 +/- 0.0295 | |
| 8 | (3,2) | L=3 'A' root 2 | 70.347388 | 69.393851 | 69.0171 +/- 0.0295 | |
| 9 | (1,3) | L=1 'A' root 3 | 72.664208 | 71.661196 | 71.1896 +/- 0.0295 | |
| 10 | (5,1) | L=5 'A' root 1 | 77.189949 | 76.242851 | 75.9521 +/- 0.0295 | |
| 11 | (2,3) | L=2 'S' root 3 | 88.206149 | 87.225441 | 86.7875 +/- 0.0295 | |
| 12 | (7,0) | L=7 'A' root 0 | 88.715954 | 87.748556 | 87.5223 +/- 0.0295 | |
| 13 | (4,2) | L=4 'S' root 2 | 89.656384 | 88.708986 | 88.3571 +/- 0.0295 | |
| 14 | (6,1) | L=6 'S' root 1 | 98.647176 | 97.698595 | 97.4193 +/- 0.0295 | |

Not targets: no J = 0 level (v = 5 lies above M = 105, so the critical_masses_v4 method is not needed here), no
L >= 8 level (first L = 8 BO threshold > 105), no 2p sigma_u level (first at 330.6).

Order of appearance (BO and BO+ad identical): (1,0) < (0,1) < (2,0) < (1,1) < (3,0) < (2,1) < (0,2) < (4,0) < (1,2)
< (3,1) < (5,0) < (2,2) < (4,1) < (0,3) < (6,0) < (3,2) < (1,3) < (5,1) < (2,3) < (7,0) < (4,2) < (6,1) < (0,4).
Restricted to S, P, D it is exactly Frolov's 1992 II list. Close pairs where the exact order could differ from BO:
(2,3)/(7,0) (BO+ad 0.52 apart, predictions 0.73 apart) and (7,0)/(4,2) (0.96 apart).

## Method

1. **BO map** (`bo_map.py`, done, 991 s): zero-energy coupling problem of `paper_data/bo_offset/coupling.py`
   (lam = M + 1/2) for L = 0..10, BO and BO+ad, all roots up to M = 105. BO roots agree with `bo_h2plus`
   count_levels bisection to < 1e-9 for six levels; BO and BO+ad level counts at M = 100 agree for J = 0..10.
2. **Frozen predictions** (`predict_map.py --freeze`, 2026-11-16 09:35, sha256
   `6e494c54232cbadcc3423309babc48fbc8fd7a1a83a5cdb2f3fff8ef99e02fba`, `PREDICTION_MAP.md`): the parameter-free
   fold of the threshold weights with m_v(R), m_r(R), non-perturbative, as in the two earlier blind tests. Band
   +/- 0.0295 for every target = largest |exact - pred| of the nine known levels. Do NOT rerun with `--freeze`.
3. **Exact** (`run_map.py`): generalization of `paper_data/thresholds_Lv/run_exact.py` to any (L, k):
   * fixed bases (mass-independent matrices), Brent root of E_k(M) - E_th(M): each crossing is a rigorous upper bound
     to M_c up to rounding; rounding-safe bound with a 1e-11 Eh margin; tol 1e-10 sensitivity;
   * compact total-degree blocks: k = 0: root-0 optimum + 0.6 alpha_X copy; k = 1: root-0 and root-1 optima;
     k >= 2: root-0, root-k optima and their geometric mean. Scales optimized at M_opt = BO+ad (x 1.03 for k >= 2), N_opt
     = 10 (L <= 2) or 8 (L >= 3);
   * nested diffuse sets (atom x X-atom relative motion, plus mirror): 5 for L = 1, 3 for L = 2, 3, 2 for L >= 4;
   * N_D chosen automatically as the largest even N with the largest check basis (N_D + 2, all diffuse sets)
     <= 8000 functions (~6 GB); compact sequence N_D - 4 .. N_D + 2; geometric extrapolation in N;
   * root-order check: the gaps of roots k - 1 and k + 1 at every M* are stored (root k - 1 must be bound);
   * binding energies at M*(1 + x), x = 5e-4 .. 0.2 (linear threshold law for L > 0).
4. **Analysis** (`analyze_map.py` -> `RESULT.md`): per-target tables, bounds, extrapolation; blind-test table
   (sha256 checked); Frolov comparison (his "bound at" mass is a rigorous upper bound our value must respect);
   complete map and order of appearance (exact vs BO vs BO+ad vs Frolov); n(M) exact vs BO vs BO+ad; muonic check
   (pp mu 8.880, dd mu 17.752, tt mu 26.585).

**L = 5, 6, 7** need new sympy term lists. `solvers/` is not modified: `tb_local.py` redirects the angular term cache
to `_cache/angular_terms_map.pkl` (seeded with a copy of the solver cache, L <= 4), and `gen_terms.py L` generates the
lists of one L into `_cache/terms_L<L>.pkl` and merges them into the main local cache under a lock (single core,
resumable from the partial `terms_L<L>.pkl`, per-L lock so parallel workers never duplicate the work). Measured: L = 5
4086 s, L = 6 9643 s (machine loaded; both merged into `_cache/angular_terms_map.pkl`); L = 7 started 09:55 and was
still running at 12:50 (terms_L7.pkl 0.46 MB, growing).

Dependence of the L >= 5 targets (J5v0, J6v0, J5v1, J7v0, J6v1) on these caches:
* J5v0, J5v1, J6v0, J6v1: need L = 5, 6, already in `_cache/angular_terms_map.pkl` (validated, below).
* J7v0: needs L = 7. If `_cache/angular_terms_map.pkl` does not yet hold them, the launcher calls
  `gen_terms.py 6` (or 7) before the first such target; this costs hours on one core (estimate 3-8 h per L, not
  measured to completion) but only blocks that worker. To move the work to the Mac Studio, copy (or commit) the whole
  `paper_data/threshold_map/_cache/` folder, including the partial `terms_L6.pkl`, `terms_L7.pkl`: generation resumes
  from them. On the Studio they can also be generated ahead of the run, in parallel and on idle cores:
  `python3 paper_data/threshold_map/gen_terms.py 6 &` and `... gen_terms.py 7 &`.
* After generation the launcher runs `validate_L.py L` once (H2+ check below, output `validate_L<L>.txt`); if it does
  not PASS, the targets of that L are skipped and reported in the log.

**Validation of L >= 5** (`validate_L.py`): the general-L tensor code was validated for L = 3, 4 in thresholds_Lv only.
Check: delta(J) = E_exact(v=0, J) - E_BO+ad(v=0, J) of H2+ (BO+ad from bo_h2plus) with references for J = 0..4 (Tang et
al. 2014 for J <= 3, our validated L = 4 value) is almost exactly linear in J(J+1) (quadratic fit residual 2e-11);
extrapolated to J = 5 it gives -6.793e-7 Eh, the new L = 5 term lists give -6.245e-7 (N = 20): difference 5.5e-8,
against a rotational energy of 3.9e-3 Eh. PASS. Honest caveat: the N = 16 -> 20 energy change at J = 5 is 1.2e-6, so the
basis error is not below ~1e-7 and the 5.5e-8 agreement is partly fortuitous; the check excludes errors larger than a
few 1e-6 (relative ~1e-3 of the rotational energy) and makes smaller ones unlikely. The S+A = unsymmetrized check of
thresholds_Lv is not repeated for L >= 5 (it needs term lists for every component pair, another ~1 h of sympy per L).
L = 6 (`validate_L6.txt`): PASS only formally at tol 1e-12 (delta +9.4e-7 vs -7.4e-7 expected; N = 24 gives no
improvement). Diagnosis: the canonical orthogonalization at tol 1e-12 drops many functions at high L (kept 2438 of
3344); at tol 1e-14 the same N = 20 basis gives delta = -7.00e-7, within 4e-8 of the extrapolation, so the L = 6 term
lists are correct. But overlap noise also grows with L: for the L = 5 H2+ basis lam_min/lam_max = -2.3e-12, and tol 1e-14
then yields a spurious root (E = -2.41). Consequences for production: (i) tol stays 1e-12 (tol 1e-10 is the sensitivity
check), so convergence in N is slower at high L and the crossings converge from above more slowly (still rigorous
bounds); (ii) `run_map.py` now logs lam_min/lam_max of every basis and raises on a spurious root (E_0 - E_th < -0.15),
so a noise root can never produce a fake crossing.

## Short tests (after the freeze; NOT results; machine loaded by another run, load 10-70)

| test | basis | size (kept) | build | project + 1 solve | max RSS | observation |
|---|---|---|---|---|---|---|
| J2v1 scales (N_opt 12) | 2 x 70 evals | ~1400 | | 2192 s total | 0.43 GB | root 1 bound at BO+ad by 2.3e-4 |
| J2v1 | N=16 + 3 diffuse | 6216 (4061) | 108 s | 138 s | 1.94 GB | at M = 27.75: E_1 - E_th = -3.2e-5, E_2 - E_th = +1.7e-6 |
| J1v2 scales (N_opt 10) | 79 + 78 evals | ~700 | | 59 s | | root 2 bound at M_opt by 1.4e-4 |
| J1v2 | N=16 + 5 diffuse | 6579 (4273) | 66 s | 112 s | 1.59 GB | at M = 39.45: E_2 - E_th = -1.5e-5 |
| J5v0 scales (N_opt 8) | 189 evals | 495 | | 598 s | | L = 5 build ~3 s for 495 functions |
| J5v0 | N=12 + 2 diffuse | 5430 (3855) | 211 s | 91 s | 1.77 GB | at BO+ad 49.35: E_0 - E_th = -2.1e-4 |

A mini end-to-end pipeline (scratchpad copy, J2v1, sizes <= 1700, all 5 stages + analyze_map.py, rerun skips
everything) ran in 58 s. Root ordering behaves as expected in all tests (roots below k bound, root k + 1 above threshold). Both of Frolov's
"bound at" masses tested (27.75, 39.45) are confirmed bound. Subset eigensolver (8 threads, loaded): 0.5 s at n = 2000,
3.5 s at 4000, 8.5 s at 5500.

## Cost (estimate from the tests; uncertain by about +/- 50 %)

A crossing on a production basis (6000-8000 functions) costs build + projection (~1-2 min) + ~15 solves (~5 s each),
i.e. ~3-5 min for L <= 2, ~6-8 min for L = 3, 4 (builds ~3x dearer per n^2), ~8-12 min for L >= 5 (L = 5 build 211 s
at n = 5430; L = 6, 7 not measured, assumed 1.5-2x L = 5). Per target (4 compact + up to 5 diffuse crossings, 5 check
crossings, 9 binding solves, scale optimization):

| targets | per target | subtotal |
|---|---|---|
| L <= 2: J1v2 J2v1 J2v2 J1v3 J2v3 | ~1 h | ~5 h |
| L = 3, 4: J3v1 J4v1 J3v2 J4v2 | ~1.25 h | ~5 h |
| L >= 5: J5v0 J6v0 J5v1 J7v0 J6v1 | ~2 h (+ scale optimization 10-30 min) | ~10 h |
| **total sequential** | | **~20 h** (plus L = 6, 7 term generation if not done: hours, one core) |

Peak memory per process: the largest basis (<= 8000 functions) needs ~6 GB at 12 n^2 doubles; measured RSS was
~6.2 n^2 doubles (1.6-1.9 GB at n = 5400-6600), so expect <= 3.5 GB per process in practice.

**Recommendation for the Mac Studio (M1 Ultra, 20 cores, 128 GB):** `WORKERS=4 THREADS=5` (memory <= ~25 GB) gives
roughly 6-8 h wall time, Frolov's three targets finishing within the first ~1.5 h. If paper_data_run6 runs at the same
time, use `WORKERS=2 THREADS=5` for each run. Leave 1-2 cores free if L = 6, 7 terms are still being generated.
On this laptop (10 cores, 16 GB): `WORKERS=1 THREADS=8`, ~20 h; run order makes partial results useful (Frolov's
three first, then by M_c).

## Launch

From the project root:

    WORKERS=4 THREADS=5 nohup caffeinate -dims scripts/paper_data_run5.sh > runs_run5.log 2>&1 &

Environment variables: `PY` (interpreter), `THREADS`, `WORKERS`, `TARGETS` (subset, priority order; default
`run_map.py --list`). Resumable: rerun the same command. Each target writes only `data/<T>_*` (parallel-safe);
`analyze_map.py` writes `RESULT.md` at the end (it can also be run any time on partial data). If a run is killed while
a validation is running, remove the stale lock `rmdir paper_data/threshold_map/data/.validate_L<L>.lock`.
Single target: `THREEBODY_BLAS_THREADS=8 python3 paper_data/threshold_map/run_map.py --target J1v2`.
Plan of a target without computing: `run_map.py --target J1v2 --dry --nocal`.

## Caveats

* Basis sizes are capped at 8000 functions; L >= 5 uses compact N <= 14 and 2 diffuse sets, k >= 2 targets 3 compact
  sets. Extrapolation errors will be larger than for J1v1/L3A/L4S; the crossings stay rigorous upper bounds.
* J2v1 scales were optimized at N_opt = 12 (the first test), the others at 10 (L <= 2) or 8 (L >= 3).
* The prediction band (+/- 0.0295) is the largest calibration residual, applied uniformly; all six L > 0 calibration
  residuals were negative (`FREEZE_NOTES.md`), so systematic failures, if any, are expected on the low side for low L, high v.
* The muonic bound-level sets in `analyze_map.py` are quoted from memory (standard facts) and should be checked against
  a review before quoting.
