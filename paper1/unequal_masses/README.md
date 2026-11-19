# unequal_masses: a blind test of the BO-offset theory with unequal nuclear masses {A+ (M), B+ (rM), e-}

*Last updated: 2026-11-16 12:21*

**STATUS: DONE.** Results: `RESULT.md`. Theory: `THEORY.md`. Frozen prediction: `PREDICTION_UNEQUAL.md`; what was seen
before the freeze: `FREEZE_NOTES.md`.

## Why this family
In the symmetric family {X+, X+, e-} the non-adiabatic offsets M_c(BO+ad) − M_c(exact) are almost constant
(0.50–0.54 for L = 0), so a "same offset as the previous level" null model also passes (referee M4).
With unequal nuclear masses the electronic problem is still H2+, but (i) the true threshold is A⁺ + (B e),
lower than the mean-atom asymptote of the 1sσg BO+ad curve by β = (r−1)/(2(1+r)) in scaled units, and
(ii) the β∂_z coupling mixes 1sσg with 2pσu, which is degenerate with it at large R. The theory
(THEORY.md) treats this non-perturbatively (two-channel g/u model, exact couplings, exact threshold, 2×2
second-order mass matrices). It predicts offsets of several units in M that change by ~2 from one v to the next:
neither persistence nor transfer of the symmetric offsets can reproduce that.

## Files
| file | content |
|---|---|
| THEORY.md | derivation (λ map, DBOC with β² t, g/u two-channel model, 2×2 mass matrices, null models) |
| ../../solvers/bo_h2plus/unequal.py | new module: masses, λ(M), thresholds, g/u coupling functions (tests: `solvers/bo_h2plus/tests/test_unequal.py`) |
| massfunc_gu.py, build_tables.py, massfunc_gu.csv | 2×2 mass matrices (Sternheimer, Σg/Σu/Πg/Πu) and table builder |
| twochannel.py | zero-energy coupling problems: BO, BOad, BOad_thr, BOad+m, 2ch, 2ch+m (+ variants) |
| predict.py | model critical masses, frozen prediction (`--freeze`), `bo_brackets.json` (brackets for run_exact), `model_roots.json` (cache) |
| PREDICTION_UNEQUAL.md, prediction_unequal.json, prediction_unequal.sha256 | FROZEN (chmod 444) |
| validate_hdplus.py, validation_hdplus.txt | HD+ v = 0, 1 at physical masses vs Bubin–Adamowicz and Baskerville 2016 |
| run_exact.py | exact runner (resumable; `--dry`, `--test N M`, `--stage`, `--list`) |
| analyze.py | RESULT.md: bounds, extrapolation, threshold-law fits, blind comparison with theory and null models |
| test_unequal_masses.py | tests (r = 1 limits, sum rules, numerics, LeanBasis vs CoupledBasis) |
| ../../scripts/paper1_run6.sh | launcher (log `runs_run6.log`, git-ignored) |

## Commands (project root)
```bash
PY=python3
$PY -m pytest -p no:cacheprovider paper1/unequal_masses/test_unequal_masses.py -q      # ~40 s
(cd solvers && $PY -m pytest -p no:cacheprovider bo_h2plus/tests -q)                   # ~7 s
$PY paper1/unequal_masses/run_exact.py --target J0v1 --dry                              # estimate
nohup caffeinate -dims scripts/paper1_run6.sh > runs_run6.log 2>&1 &                     # PRODUCTION (user)
tail -f runs_run6.log
$PY paper1/unequal_masses/analyze.py                                                    # RESULT.md
```
Stop: `pkill -f paper1_run6.sh; pkill -f unequal_masses/run_exact.py`. Resume = same command.
Do NOT rerun `predict.py --freeze` (it refuses anyway while the sha256 file exists).

## Frozen predictions (r = 2, light-nucleus mass M; 09:48, sha256 7ac806bf...)

| level | BO | BO+ad (own asymptote) | BO+ad (exact threshold) | 2ch (no masses) | **theory 2ch+m** | null NS | null NT |
|---|---|---|---|---|---|---|---|
| J0v1 | 8.3528 | 7.5496 | 14.8883 | 11.9591 | **11.713 ± 0.030** | 7.146 | 14.485 |
| J0v2 | 23.4175 | 22.6099 | 34.1938 | 29.3779 | **29.089 ± 0.030** | 22.230 | 33.814 |
| J0v3 | 45.6767 | 44.8621 | 61.2101 | 53.7825 | **53.474 ± 0.030** | 44.489 | 60.837 |
| J0v4 | 74.8897 | 74.0691 | 95.8874 | 84.7023 | **84.386 ± 0.030** | 73.699 | 95.518 |
| J1v0 | 4.1948 | 3.4264 | 6.5121 | 5.7314 | **5.582 ± 0.037** | 3.155 | 6.240 |
| J2v0 | 9.5965 | 8.8330 | 11.7526 | 11.3059 | **11.148 ± 0.030** | 8.616 | 11.536 |

r = 1.5 (secondary): J0v1 11.812, J0v2 30.610, J0v3 57.310, J0v4 91.429, J1v0 5.408, J2v0 11.453 (PREDICTION_UNEQUAL.md).
Lenient criterion (also frozen): |exact − theory| ≤ 0.10. Theory − null differences: NS 2.4–10.7, NT 0.39–11.1;
predicted offsets from BO+ad change by −2.32, −2.13, −1.70 between consecutive v (persistence predicts 0).

## Earlier estimates (BF92)
Added 2026-11-16 after reading Bishop and Frolov, J. Phys. B 25, 3049 (1992) in full. Their Tables 3-8 give low-precision thresholds of the lowest P and D levels for unequal nuclear masses (700 or 800
exponential functions, no error bars; each threshold is an estimate between a tabulated bound and a tabulated unbound point). With
m_Z = m_e and m_X = 2 m_Y = 2 lambda (Tables 3 and 4; our r = 2 with the light-nucleus mass M = lambda, mapping assumed, to be
checked against their definitions):

| target | BF92 estimate | last bound / first unbound in their table | frozen theory (r = 2) |
|---|---|---|---|
| J1v0 (lowest P) | 5.425 (Table 3) | 5.5 / 5.4 | 5.582 +- 0.037 |
| J2v0 (lowest D) | 11.20 (Table 4) | 11.5 / 11.0 | 11.148 +- 0.030 |

For r = 3 (m_X = 3 m_Y; not a frozen target): P 5.475 (bound 5.5, unbound 5.45), D 10.65 (bound 10.75, unbound 10.5); for m_X infinite:
P 5.30, D 9.375. **To be compared after the run; no prediction file was changed.** Flag: a level found bound at M = 5.5 by a
variational calculation implies an exact threshold below 5.5, which would lie 0.08 below the frozen J1v0 theory value (2.2 of its
quoted +- 0.037); if the mapping above is right and the run confirms a J1v0 threshold above 5.5, the BF92 table point would have to be
explained (different state or sector, mass convention, or a misprint: the notes record several misprints in their tables). The D estimate
is consistent with the frozen value. Do not draw conclusions before the exact run.

## Validation (before the freeze; no crossing search)
* HD+ at physical masses (`validation_hdplus.txt`, unsymmetrized basis, weights (1,1,2), N = 36, size 4750):
  E(v=0) = −0.5978979686093 vs Baskerville 2016 −0.59789796860 (−9e−12); vs Bubin–Adamowicz 1500 CECG
  −0.597897968588 (ours 2e−11 lower, variational); E(v=1) = −0.5891818294756 vs −0.589181829513 (+3.7e−11,
  still converging). 274 s, 1.7 GB.
* r = 1 limits of all model codes vs paper1/bo_offset (1e−10) and massfunc (0.5281); unsymmetrized basis at r = 1
  equals the symmetric 'S' energies (M = 10, −0.45404517). An exact r = 1 crossing (target R1J0v1, expected
  9.525634) is part of the production run.

## Timing and memory (measured 2026-11-16 under a machine load of 30–80 from another job: inflated)
* Stage 1 (two Nelder–Mead scale optimizations, N = 12): 895 s measured (unloaded estimate 1–3 min).
* `--dry J0v1`: basis 4095: build+project 208 s, one eigensolve 75 s, peak RSS 0.90 GB (= 6.7 n² doubles).
  The printed total (16 h for J0v1) uses these loaded numbers. Unloaded reference (critical_masses_v4,
  symmetric n = 3330): build + one solve 27 s. Scaled from that: ~3 min per crossing at n ≈ 5000, ~15–20 min
  for the largest check basis; **~0.7–1.2 h per L = 0 target, ~0.5–1 h per L > 0 target, total ~6–9 h**
  for all 8 targets on an idle machine. Run `--dry` again on an idle machine for a real estimate.
* Largest bases: L = 0 diffuse-enlarged check 9543 (J0v4: 10512) → peak ~5–6 GB; L = 1/2 ≤ 8000 → ≤ 3.5 GB.
  All below the 8 GB limit. Do not run other heavy jobs in parallel.
