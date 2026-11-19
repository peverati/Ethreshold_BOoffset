# massfunc: first-principles nonadiabatic mass functions m_v(R), m_r(R) for {X+, X+, e-}

*Last updated: 2026-11-14 17:08*

**Status: fine run done 2026-11-14 (finished 16:59), `results_fine.csv`, 2597 R points; see `analysis.txt`, `fig_massfunc.png`, `fig_fold.png`.** Convergence 3e-11 (m_v), 1e-12 (m_r). m_v(R_e) = 0.346552, m_r(R_e) = 0.075422. Parameter-free non-adiabatic offsets (non-perturbative fold): 0.5281, 0.5124, 0.5017 (L=0 v=1-3; exact 0.5375, 0.5062, 0.4970), 0.3328 (L=1; exact 0.3623), 0.2778 (L=2; exact 0.2892).

Written 2026-11-14. Computes the R-dependent vibrational mass m_v(R) and rotational mass m_r(R) of H2+-like {X+, X+, e-} from second-order nonadiabatic perturbation theory. It then uses them to predict, **with no fitted parameter**, the nonadiabatic part of the BO+ad critical-mass offsets studied in `paper1/bo_offset/README.md`.

* Theory: `THEORY.md`. The project convention is λ = M + 1/2, r = μ_e R, and M → M + m. In it, m_v = 4⟨∂_Rφ0|Q(h−E0)^{-1}Q|∂_Rφ0⟩ (Σg only) and m_r = 4⟨L_yφ0/R|(h−E0)^{-1}|L_yφ0/R⟩ (Πg only). The separated-atom limit m_v(∞) = m_r(∞) = 1/2 is derived there.
* Method: the Sternheimer/Hylleraas problem is solved in a prolate-spheroidal Laguerre × Legendre basis. All H and S integrals are done by exact Gauss quadrature. The right-hand sides come from the exact bo_h2plus wavefunction.

## Commands (run from the project root)

```bash
PY=python3
# tests (~2 s)
$PY -m pytest -p no:cacheprovider paper1/bo_offset/massfunc/test_massfunc.py -q
# plan and time estimate only
$PY paper1/bo_offset/massfunc/run.py --dry            # default grid
$PY paper1/bo_offset/massfunc/run.py --dry --fine     # fine grid
# production (unattended, logs to runs_massfunc.log at the project root)
nohup caffeinate -dims scripts/massfunc_run.sh > runs_massfunc.log 2>&1 &
#   or the longer, higher-accuracy version:
nohup caffeinate -dims scripts/massfunc_run.sh --fine > runs_massfunc.log 2>&1 &
# follow it
tail -f runs_massfunc.log
# stop
pkill -f massfunc_run.sh; pkill -f massfunc/run.py
```

* **Resume:** rerun the same launcher command. R values already in `results.csv` (or `results_fine.csv`) are skipped, because one row is written and flushed per R.
* **Analysis only** (after a run, or on a partial file): `$PY paper1/bo_offset/massfunc/analyze.py [--results paper1/bo_offset/massfunc/results.csv] [--no-nonpert]`.
* **Single points:** `$PY paper1/bo_offset/massfunc/run.py --R 2 5 10 --out /some/scratch.csv` or `$PY paper1/bo_offset/massfunc/massfunc.py 2.0 6.0`.

## Expected cost (this Mac, 10 cores, 4 BLAS threads, machine load about 8 during the tests)

| run | R points | basis levels per R | wall time | peak memory |
|---|---|---|---|---|
| default | 650 (0.02–60 bohr; step 0.02 up to R = 10) | 0 and 1 (value = level 1) | about 8–10 min | < 1 GB |
| `--fine` | 2597 (4× denser) | 1 and 2 (value = level 2) | about 3 h (`--dry` estimate 172 min) | about 1.5 GB |

The calculation is much cheaper than anticipated (0.1–2 s per R). The default run is therefore shorter than 30 min. The test results below show that it is already converged to ≤ 1e−10. `--fine` mostly buys a denser grid and larger-basis cross-checks; it does not change the physics. Threads: `MASSFUNC_THREADS` (default 4).

## Outputs

* `results.csv` / `results_fine.csv`, one row per R:
  * R (unit-mass bohr, = μ_e R_phys);
  * m_v, m_r;
  * dm_v, dm_r: |difference between the two basis levels|;
  * the lower-level values and A_par, A_perp;
  * E0_exact and dE0_basis (basis ground-state error);
  * fpar_ratio, fperp_ratio: basis completeness for the right-hand sides vs bo_h2plus D(R), filled every 25 points;
  * lowest excited Σg and lowest Πg energies, basis sizes, seconds, time stamp.
* After `analyze.py`:
  * `analysis.txt`;
  * `fold_results.csv` (predicted vs exact offsets);
  * `threshold_weights.csv`: the zero-energy weight densities w_rad(r), w_cen(r) of the BO+ad threshold states, exported from the parent folder's `coupling.py` (offset = ∫m_v w_rad + ∫m_r w_cen);
  * `fig_massfunc.png` (m_v, m_r and the L = 0 weights);
  * `fig_fold.png` (predicted vs exact offsets).

## Test results (2026-11-14, short tests only; the production run was NOT done)

* **Basis convergence:** at R = 0.5, 2, 6 and 15, bases from (n_ξ, n_η) = (16, 10+1.3R) up to (40, 30+1.3R), together with exponent changes ×0.7 and ×1.4, agree to ≤ 1e−12. Small R (0.02–0.1) converges to ~1e−10. |E0_basis − E0_exact| ≤ 1e−12 for R ≥ 0.5 and ≤ 2e−9 at R = 0.05.
* **Values:**

| R | m_v | m_r |
|---|---|---|
| 0.5 | 0.250241 | 0.008336 |
| 1.0 | 0.303427 | 0.025804 |
| R_e = 1.99719 | **0.346552** | **0.075422** |
| 3 | 0.425376 | 0.143557 |
| 4 | 0.537161 | 0.226210 |
| 6 | 0.660321 (maximum near here) | 0.382848 |
| 8 | 0.591693 | 0.461617 |
| 10 | 0.535736 | 0.486523 |
| 20 | 0.502279 | 0.498871 |
| 40 | 0.500282 | 0.499859 |

* **Data comparison:** m_v(R_e) = 0.34655 vs the κ⁶ value 0.3465 ± 0.0003 (0.2σ). m_r(R_e) = 0.07542 vs the J-spacing value 0.0754 ± 0.001.
* **Large R:** m − 1/2 ∝ R^{−3} (ratio between R = 20 and R = 40 is 8.1). This confirms the 1/2 limit.
* **Hump:** m_v(R) overshoots 1/2, with a maximum of about 0.66 near R ≈ 6 bohr. This confirms the inference in the parent README. m_r rises monotonically.
* **Preliminary fold** (coarse 67-point test grid, not the production grid). Predicted nonadiabatic BO+ad offsets, first order:

| state | predicted | exact |
|---|---|---|
| L = 0, v = 1 | 0.526 | 0.5375 |
| L = 0, v = 2 | 0.512 | 0.5062 |
| L = 0, v = 3 | 0.502 | 0.4970 |
| L = 1 | 0.331 | 0.3623 |
| L = 2 | 0.278 | 0.2892 |

  The nonperturbative q-variant agrees within 0.003. The residuals are ≤ 0.012 for L = 0 and L = 2, and 0.031 for L = 1. At M_c = 4.2 the neglected O(1/M) terms are expected to matter most there.

## Files

* `massfunc.py`: core (basis, matrices, Sternheimer, `mass_functions(R)`).
* `run.py`: driver.
* `analyze.py`: analysis.
* `test_massfunc.py`: tests.
* `THEORY.md`: derivation.
* Launcher: `scripts/massfunc_run.sh`.
