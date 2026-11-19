# M = 1e4 convergence (L = 0, 1, 2)
*Last updated: 2026-11-15 15:22*

**Why:** at M = 1e4 the mass-scan levels are the least converged points of Paper 1. Their current errors are L = 0 2.2e-9, L = 1 2.1e-8 and L = 2 1.0e-7; the L = 2 error leaves d2/d1 = 3.0009 ± 0.002 unresolved (see `solvers/threebody/README.md`, rotational table). The target is about 1e-10 for L = 0 and better than 1e-9 for L = 1, 2.

**What:** `run.py` continues the nested fixed-scale sequences of `paper_data/mass_scans` with the same alphas, which it reads from `exact_L0.csv` and `exact_rot.csv`. The new basis sizes are:

| Level | Tier 1 | Tier 2 |
|---|---|---|
| L = 0 | N = 52, 56, 60 | N = 64 |
| L = 1 | N = 42, 46 | N = 50 |
| L = 2 | N = 38, 42 | none |

Tier 2 has the largest matrices (12–13k) and may swap on 16 GB.

**Checkpointing:** one row per (level, N) is written to `results.csv` as soon as it finishes, and the log goes to `run_log.txt`. Rerunning skips the rows already done, so a crash loses at most the N in progress. Each level's `SUMMARY` line in the log gives the extrapolation over the old and new N together.

**Run:** via `scripts/paper_data_day_run.sh`, which also rebuilds `solvers/threebody/results/codata2012/`. `python run.py --dry` prints the plan without computing anything.

**Status: done 2026-11-14 08:29–09:30 (1 h total, 2–11 min per N, no swapping).**

| Level | New N (size) | E at largest N | Extrapolated | Error | Old extrapolation |
|---|---|---|---|---|---|
| L = 0 | 52–64 (12529) | −0.60032852464703 | −0.60032852464759 | 8.3e−13 | −0.6003285245581 ± 2.2e−9 |
| L = 1 | 42–50 (12051) | −0.60027897259580 | −0.60027897267350 | 1.4e−10 | −0.6002789831763 ± 2.1e−8 |
| L = 2 | 38–42 (11132) | −0.60017993272494 | −0.60017994634621 | 1.4e−8 | −0.6001798551815 ± 1.0e−7 |

The old values were within their error bars. The new spacings are d₁ = 4.955197e−5 ± 1.4e−10 and d₂ = 1.485783e−4 ± 1.4e−8, so **d₂/d₁ = 2.9984 ± 0.0003**, against the BO value 2.9985. The anomalous 3.0009 is gone: the ratio sits below the rigid-rotor 3 like every other mass. M·d₁ = 0.49552, against 2⟨1/R²⟩ = 0.4956. E − E_BO+ad = −4.005e−8, with its error now set by the 10-digit BO+ad value. L = 2 is still the weakest level; N = 46 (size about 14.5k) would be the next step.

## Run 3: L = 2, N = 46 (done 2026-11-15 14:16–14:33, 1022 s, lean path)

**Result:** E(N=46) = −0.60017993973803; extrapolation over N = 26–46: **E(L=2) = −0.60017994134773 ± 5.0e−9** (ρ = 0.187), within the previous −0.60017994634621 ± 1.4e−8. Hence d₂ = 1.485833e−4 ± 5e−9 and **d₂/d₁ = 2.99853 ± 0.00010**, equal to BO (2.99853). `paper_data/mass_scans/analyze.py` picks the new row up automatically (`spacings.csv` updated; `fits.txt` unchanged).

### Setup notes (as planned)

* **L = 2, N = 46, size 14400** (`JOBS` entry `("L2", 2, "S", [46])`). Rows already in `results.csv` are skipped (checked with `--dry`: only L2 N=46 is TO DO).
* **Memory.** The standard path (`CoupledBasis` keeps 5 raw operator matrices, the projected copies and LAPACK work) was measured at about 13 n² doubles of peak footprint (2.03 GB at n = 4352). At n = 14400 (1.66 GB per dense matrix) that is about 22 GB: **too risky on 16 GB** (it would swap heavily).
* **Safer option, used:** `solve_lean()` in `run.py` (pairs listed in `LEAN`). It diagonalizes S alone, assembles H as a single matrix, projects and diagonalizes. Same function space, tol and energy: N = 26 and 30 agree with the standard path to 6e-14 and 2e-13. Measured peak footprint at N = 30: 0.53 GB (standard 2.03 GB), i.e. about 2.5–4 n², so **about 4–7 GB at N = 46**.
* **Expected time:** 20–40 min with 8 BLAS threads (N = 42 took 679 s on the standard path; cost scales roughly as n³, factor 2.2).
* **Launch:** part of `scripts/paper_data_run3.sh` (runs first). Alone: `THREEBODY_BLAS_THREADS=8 python3 run.py`.
* If even the lean path is too heavy on a loaded machine, delete the `[46]` job line; nothing else depends on it.
