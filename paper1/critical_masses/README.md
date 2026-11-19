# Critical masses of the symmetric S states of {X+, X+, e-}

*Last updated: 2026-11-16 08:18*

Written 2026-11-12. The system is two identical particles of mass M and charge +1 plus one particle of mass 1 and charge −1 (atomic units), with L = 0 and the spatial wave function symmetric under exchange of the two X. M_ck is the mass at which the (k+1)-th symmetric S state (root k; root 0 is the ground state) becomes bound below E_th(M) = −M/(2(M+1)).

## Results

| k | M_ck (non-BO, extrapolated) | upper bound, margin 0 | rounding-safe upper bound (margin 1e−11 Eh) | BO | BO + adiabatic | BO − exact | (BO+ad) − exact |
|---|---|---|---|---|---|---|---|
| 1 | **9.525634 ± 0.000004** | 9.525652 | 9.52586 | 11.0814566 | 10.0631790 | 1.556 | 0.538 |
| 2 | **29.63939 ± 0.00001** | 29.639433 | 29.64016 | 31.1677743 | 30.1455552 | 1.528 | 0.506 |
| 3 | **59.31866 ± 0.00002** | 59.318722 | 59.32056 | 60.8466744 | 59.8156817 | 1.528 | 0.497 |

* **Rigorous.** For M at or above the "rounding-safe" bound, a variational root lies at least 1e−11 Eh below E_th. Observed rounding is about 1e−13 Eh, so this proves that the state is bound there. The same statement with margin 0 (the "upper bound, margin 0" column) is rigorous only up to rounding. It needs a margin because the fixed-basis root crosses threshold with a very small slope (1.0e−8 Eh per mass unit for k=1, 1.9e−9 for k=2, 5e−10 to 8e−10 for k=3). Every M_ck is therefore proven to be below the BO and BO+adiabatic predictions.
* **Extrapolated, not rigorous.** The lower side, and so the central values, come from the threshold-law fits below. A variational calculation cannot prove that a state is unbound.
* The BO and BO+adiabatic values come from `solvers/bo_h2plus`, unchanged; see `bo_critical_masses.csv`.
  * BO+ad is counted at that model's own asymptote, −μ_e/2 + μ_e²/(4M).
  * Requiring BO+ad levels to lie below the *exact* threshold instead gives 12.136, 32.491, 62.630. This is an inconsistent criterion, because the BO+ad asymptote lies about 1e−3 Eh above E_th at M ≈ 10. It is listed only for completeness.
  * BO+ad has its first level bound for every M ≥ 0.05 tested, so it has no first critical mass in that range.
* **Observation (not explained).** The exact critical masses lie below BO by an almost constant ΔM ≈ 1.53–1.56, and below BO+ad by ≈ 0.50–0.54, for all three states.

## Threshold law

**Expected law.** Near threshold the state is X+ weakly bound, in an s wave, to the neutral atom (X e)1s. At long range the ion-atom interaction is the polarization tail −C4/R⁴. The atom is in a nondegenerate 1s state, so there is no 1/R² dipole term, and the g/u charge-exchange splitting decays exponentially. This is short-range scattering in the threshold sense: the scattering length diverges linearly in the coupling, so κ = sqrt(2 μ E_b) is linear in x = M − M_c to leading order, and **E_b ∝ (M − M_c)²**. Here μ is the ion-atom reduced mass M(M+1)/(2M+1). For R⁻⁴ the effective-range expansion has non-analytic κ³ ln κ-type terms, which enter only at higher order. A Coulomb-tail system (critical nuclear charge of H⁻-like ions) would instead give a linear law, E_b ∝ x.

**Test (data in `fits_state*.txt`).** The data are 12 masses per state, from M* + 0.001 to M* + 1, with E_b from 2e−10 (state 3: 3e−12) to 2e−4 Eh.

* A free power law E_b = A x^p on the window x ≤ 0.05 gives:
  * k = 1: p = 2.0012 ± 0.0002
  * k = 2: p = 2.0044 ± 0.0008
  * k = 3: p = 2.0057 ± 0.0007

  p is 2. The small excess grows with the window width, as expected when a c2 x² term in sqrt(E_b) is absorbed into p. A linear law, p = 1, is excluded.
* sqrt(E_b) = c1 x + c2 x² + c3 x³ fits within relative rms 1e−7 (k=1), 6e−8 to 3e−7 (k=2) and 5e−6 (k=3) for x ≤ 0.1. The fitted coefficients are:

| k | c1 | c2 | c3 |
|---|---|---|---|
| 1 | 0.013737 | 4.6e−4 | −1.6e−3 |
| 2 | 0.004131 | 4.4e−4 | −4e−4 |
| 3 | 0.001666 | 2.0e−4 | −1e−4 |

**Error bars of M_c.** The quoted values are the order-3 fits on x ≤ 0.1. The error bar is the full spread of all order-2 and order-3 fits on windows x ≤ 0.05 to 0.3, plus the free-power-law fits. For k=1 these span 9.525632 to 9.525637. The extrapolation distance is only about 0.001 in M, from the lowest data point to M_c. Basis errors are negligible next to the model spread: E_b is converged to about 1e−13 Eh absolute (see below), which moves M_c by less than 1e−7.

**Cross-check.** The fixed-basis crossing masses M*(α_s) decrease monotonically as the basis becomes more diffuse:

| k | α_s = 0.01 | 3e−3 | 1e−3 | 3e−4 | 1e−4 |
|---|---|---|---|---|---|
| 1 | 9.52830 | 9.52644 | 9.52591 | 9.52572 | 9.52566 |
| 2 | 29.64449 | 29.64094 | 29.63991 | 29.63955 | 29.63944 |
| 3 | 59.32764 | 59.32139 | 59.31957 | 59.31893 | 59.31875 |

They approach the fitted M_c from above. α_s = 3e−5 is worse than α_s = 1e−4, because with n_d = 40 that block loses resolution at intermediate R. A naive linear extrapolation in α_s gives 9.52569 to 9.52572 for k=1, slightly above the fit. This finite-size extrapolation is cruder than the fit and is not used for the quoted value.

## Size of the state

As M → M_c, ⟨r_XX⟩ approaches the zero-range halo value 1/(2κ); the ratio ⟨r_XX⟩·2κ is tabulated in `fits_state*.txt` and plotted in `fig_state*.png`.

| k | E_b at the lowest data point | ⟨r_XX⟩ | ratio ⟨r_XX⟩·2κ |
|---|---|---|---|
| 1 | 2.0e−10 | 11 209 bohr | 1.0007 |
| 2 | 1.9e−11 | 20 948 bohr | 1.0009 |
| 3 | 1.2e−11 | 18 577 bohr | 1.0003 |

* The ratio grows to about 1.5 at x = 1, where finite-range effects matter.
* ⟨r_Xe⟩ ≈ ⟨r_XX⟩/2: one X sits in the atom and the other is far away.
* ⟨r_XX⟩ diverges as 1/x, and the rms value as sqrt(2) times the mean, as for exp(−κR)/R.
* The single point at E_b = 3.3e−12 for k=3 has ratio 0.984. There 1/κ ≈ 7e4 bohr and the basis is too small. That point is not converged: its ladder spread is 0.4% of E_b.

## Method

**Basis.** The new module `solvers/threebody/threebody/threshold.py` builds a union of perimetric-Laguerre blocks, each with its own three scales. Exchange symmetry is imposed as φ = F + PF, where P maps (s0, s1, s2) to (s1, s0, s2), so a block may have α0 ≠ α1. All integrals are exact:

* integrals between different scales are done with Gauss-Laguerre quadrature (the integrand is a polynomial times an exponential);
* derivatives use the Laguerre derivative matrix.

H(M) = T_X/M + T_e + V_01 − V_Xe, with mass-independent matrices, so a mass scan only re-diagonalizes. The generalized eigenproblem is solved by canonical orthogonalization at tol = 1e−12.

**Basis used for state k.** The basis is fixed and independent of M:

* **Compact blocks.** Total degree N = 16, scales (α_X, α_e):
  * ground-state optimum
  * root-k optimum (package `optimize_alphas`, N = 12)
  * for k ≥ 2, an intermediate set
* **Two diffuse blocks.** Scales (μ_a, α_s, μ_a) and (μ_a, 0.03, μ_a), with μ_a = M_ref/(M_ref+1).
  * a0 + a2 ≤ 4: the atom and its polarization, in s0 and s2
  * a1 ≤ 40 and ≤ 30 respectively: the ion-atom coordinate s1 ≈ 2R
* **Ladder.** α_s ∈ {1e−2, 3e−3, 1e−3, 3e−4, 1e−4, 3e−5}.
* **Size.** 2130 functions (k=1) and 2655 functions (k=2, 3).

At each M the lowest root over the ladder is taken. Each root is a Hylleraas-Undheim-MacDonald upper bound, so E_b is a lower bound on the exact binding.

**Convergence evidence** (`convergence_state*.csv`, `crossings_state*.csv`). Each check was made at three masses (M* + 0.01, M* + 0.1, M* + 0.5) by enlarging one part of the basis:

| change to the basis | change in E_k (all states) |
|---|---|
| compact N 16 → 20 | ≤ 1.3e−13 Eh |
| diffuse a0 + a2 ≤ 4 → 6 | ≤ 2e−13 Eh |
| n_d 40/30 → 60/45 | ≤ 1.1e−12 Eh |
| orthogonalization tol 1e−12 → 1e−10 | ≤ 6e−12 Eh (smaller space, higher energy, as expected) |

* **Ladder spread.** The best and second-best α_s agree to 1e−14 to 1e−13 Eh, except at the lowest points.
* **Earlier check at M = 10** (exploration, not saved as a table). Four unions gave the same E1 − E_th = −4.23763475e−5 to all printed digits: C16+C16+D, C20+C16+D, C20+C20+D, and C24+C20+D. The previous single-scale value was −3.8e−5 at N = 24.

**Warnings.**

* **tol = 1e−14 is unsafe.** With tol = 1e−14 a spurious root appeared at k=2 (M = 29.649), 7e−7 Eh below the converged value. This is rounding from a nearly singular overlap. It is why every value here uses tol ≥ 1e−12, and why the rigorous bound uses a margin.
* **Nesting can fail at rounding level.** Canonical orthogonalization drops near-dependent directions, so a nominally larger basis is not exactly a superset. For example, compact N=20 gave a crossing M* 1.4e−5 *higher* than N=16 for k=3, which corresponds to about 1e−13 Eh.

## Files

* `bo_critical_masses.csv`: BO and BO+adiabatic critical masses.
* `energies_state{k}.csv`: E_k − E_th for every M and α_s.
* `observables_state{k}.csv`: best E, E_b, κ, ⟨r_XX⟩, rms r_XX, ⟨r_Xe⟩, dE/dM.
* `crossings_state{k}.csv`: fixed-basis crossing masses, sizes, and timings.
* `convergence_state{k}.csv`: basis and tol checks.
* `upper_bounds.csv`: crossings with margins 0, 1e−11 and 1e−10, at tol 1e−12 and 1e−10, with the slope at M*.
* `fits_state{k}.txt`, `fig_state{k}.png`: fits and plots.
* `basis_state{k}.json`: basis definition.
* `run_log.txt`, `PROGRESS.md`: run log and working notes.

**Code** (all new files; no existing file was changed):

* `solvers/threebody/threebody/threshold.py`
* `solvers/threebody/tests/test_threshold.py`: 5 tests
  * mixed 1D integrals against quadrature
  * reproduces the package basis exactly
  * symmetrization against an unsymmetrized union
  * variational union
  * Hellmann-Feynman dE/dM
* `solvers/threebody/examples/critical_mass.py`: compute.
* `solvers/threebody/examples/critical_mass_fit.py`: fits.
* `solvers/threebody/examples/critical_mass_bounds.py`: margin bounds.
* `solvers/threebody/examples/critical_mass_bo.py`: BO and BO+ad.

**Reproduce** (from `solvers/threebody`, with `THREEBODY_BLAS_THREADS=3`):

```
python examples/critical_mass.py ../../paper1/critical_masses 1 2 3
python examples/critical_mass_fit.py ../../paper1/critical_masses
python examples/critical_mass_bounds.py ../../paper1/critical_masses
cd .. && python threebody/examples/critical_mass_bo.py ../paper1/critical_masses
```

**Run times** on a shared, at times heavily loaded 10-core Mac (load average up to about 130):

| step | time |
|---|---|
| k=1 | 940 s |
| k=2 | 1068 s |
| k=3 | 863 s |
| margin bounds | about 15 min |
| BO | 3 s |
| exploration | about 5 min |

## What remains uncertain

* **Lower side of M_c.** It is extrapolated: it rests on the analytic threshold law, which is standard for short-range s-wave binding and fits the data to 1e−7 relative. It is not proven.
* **Ladder.** The near-threshold behaviour at x < 0.001 was not sampled, and the most diffuse ladder member (3e−5) is resolution-limited.
* **Fixed compact scales.** These were optimized at a single M per state. Since E_k is converged to about 1e−13 this does not matter, but the scales were not re-optimized.
* **Only S states with J = 0 and symmetric spatial exchange.** Antisymmetric states and J > 0 states were not studied.
* **Literature (corrected 2026-11-16).** Frolov and co-workers located these thresholds earlier at lower precision: v=1 10.045, v=2 30.70, v=3 61.135 (Frolov, J. Phys. B 25, 3059 (1992)), 3 to 5.5% above ours. See docs/literature/frolov_1992_notes.md and paper1/RESULTS_SUMMARY.md (PRIORITY block).
