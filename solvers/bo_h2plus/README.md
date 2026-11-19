# bo_h2plus: exact Born-Oppenheimer reference for {X+, X+, e-}

*Last updated: 2026-11-17 20:42*

> **Mass constants (2026-11-17):** Paper 1 quotes physical-mass values at CODATA 2022 (m_p/m_e = 1836.152673426, m_μ/m_e = 206.7682827); they are in `paper_data/codata2022/` (README, results.csv) and `solvers/threebody/results/codata2022/`. Examples and tables in this README that name other masses (CODATA 2018 1836.15267343, 206.768283; CODATA 2010 1836.15267245, 206.7682843, used for the King et al. comparison and called "CODATA 2012" in older notes) are kept as documentation of those runs; they differ from the CODATA 2022 values only in the last digits.

This package computes the **Born-Oppenheimer (BO) prediction** for the three-particle Coulomb family made of two identical positive particles of mass M (in electron masses) and one electron, for any M. Examples: M = 1836.15267245 (CODATA 2012) is H2+; M of about 1 is a Ps⁻-like system, where BO is not meaningful but still well defined. All quantities are in atomic units unless a name ends in `_cm`.

Written 2026-11-12. Python 3 with numpy, scipy and threadpoolctl (pytest is needed for the tests; sympy was used only to derive the recurrence). Nothing was installed.

## Equations

Jacobi coordinates: R is the X–X vector and ρ is the electron position relative to the X–X midpoint.

    H = -(1/2mu_N) Lap_R - (1/2mu_e) Lap_rho + 1/R - 1/|rho+R/2| - 1/|rho-R/2|,
    mu_N = M/2,  mu_e = 2M/(2M+1)

The coordinates are orthogonal, so there is **no mass-polarization term**: the midpoint is the nuclear centre of mass.

* **Electronic problem** (clamped nuclei, electron mass 1): U(R) = E_el(R) + 1/R. For mass mu_e the potential is exactly V_BO(R; M) = mu_e U(mu_e R). The tests also check this numerically through the equivalence "electron mass mu ⇔ nuclear charge Z = mu, energy/mu".
* **Nuclear problem**:

      -(1/2mu_N) u'' + [V_BO(R) (+ E_ad(R)) + J(J+1)/(2 mu_N R^2)] u = E u

* **Thresholds**:
  * BO: −mu_e/2 = −M/(2M+1).
  * Exact two-body limit X + (Xe): −M/(2(M+1)).
  * The two differ by M/((2M+1)(2M+2)), about 1/(4M). For H2+ that is 1.36e-4 hartree; for M = 1 it is 0.083 hartree.
  * With the adiabatic correction the BO threshold becomes −mu_e/2 + mu_e²/(4M). This equals the exact limit to first order in 1/M only.
* **Useful exact property**: dividing H_BO by mu_e gives a unit-mass problem with effective mass m_eff = mu_N/mu_e = (2M+1)/4. So the J=0 level count n(M) depends on M only through m_eff, and E_v(M) = mu_e · ε_v(m_eff).

## Methods

1. **Electronic** (`electronic.py`) uses prolate spheroidal separation.
   * Angular equation: expansion in normalized Legendre polynomials (even l for 1sσg, odd l for 2pσu). This gives a symmetric tridiagonal eigenproblem for the separation constant λ.
   * Radial equation: Jaffé expansion X = (ξ+1)^σ e^{−pξ} Σ g_n t^n, with t = (ξ−1)/(ξ+1) and σ = R/p − 1. The three-term recurrence was **derived symbolically with sympy**:

         (n+1)² g_{n+1} + [λ − (2n² + 4np + 2n + p² + 2p + 1 − (a/p)n − a − a/2p)] g_{n} + (n−1−σ)² g_{n−1} = 0,  a = 2ZR

     Truncated, it becomes a symmetrizable tridiagonal eigenproblem.
   * The energy E = −2p²/R² comes from λ_rad(p) = λ_ang(p), solved by Brent's method.
   * One energy takes about 1–3 ms. dE/dR uses a Richardson-extrapolated central difference.
   * `equilibrium(mu)` gives R_e, E_min and k using the exact mass scaling.
   * The 2pσu state is also available (`state='2pu'`).
2. **Adiabatic (diagonal) correction** (`adiabatic.py`):

       E_ad(R; M) = mu_e² D(mu_e R)/M,   D(R) = <∂_R ψ|∂_R ψ>_ρ + <L_x² + L_y²>/R²

   * The derivative at fixed ρ is built from the analytic chain rule in spheroidal coordinates.
   * The (ξ, η) quadrature uses ξ = cosh a, η = cos b with Gauss-Legendre nodes.
   * It is an option (`adiabatic=True`) and is never included unless requested.
3. **Tables** (`potential.py`): U(R) and D(R) are interpolated with piecewise Chebyshev polynomials on 15 geometric segments over [0.01, 200] bohr.
   * The tables are cached in `_cache/*.npz`. You can delete these files; they are rebuilt in about 1 s (E_el) and 25 s (D).
   * Beyond 200 bohr: −1/2 − 9/(4R⁴) − 15/(2R⁶).
   * Below 0.01 bohr: −2 + c R² (that region is irrelevant there because 1/R dominates).
4. **Nuclear** (`nuclear.py`):
   * Grid: sinc-DVR (the Colbert–Miller radial form, with φ(0)=0) in the mapped variable x, where R = R_s sinh x. The map is uniform near R = 0 and logarithmic at large R. This yields a symmetric problem with a diagonal metric.
   * Step: dx = π/(ppw · k_max), where k_max is the largest local wavenumber in x. It therefore adapts to M: about 100 points for M ~ 1–1000, about 600 for M = 1e4, about 2000 for M = 1e5.
   * Box: R_max is set by the length β = sqrt(2 mu_N C4) of the −C4/R⁴ tail (C4 = 9/(4 mu_e³)). It is then **enlarged automatically** until two conditions hold: the DVR finds as many levels as the independent node count, and R_max ≥ 25/κ for the least-bound level.
   * Each call solves again on a second grid (ppw/1.4 and a box twice as large) and stores the difference per level as `err_E` and `err_R`.
   * **Level count** `count_levels`: the number of nodes of the threshold-energy solution, computed independently of the DVR. It uses Numerov on the same map, started from a 12-term Frobenius series at R = 0, out to R_b ≫ β. Beyond R_b it continues analytically: for J=0 with the exact zero-energy −C4/R⁴ solution u = R sin(β/R + φ), and for J>0 with a R^{J+1} + b R^{−J}.
   * Outputs per level: E, binding, ⟨R⟩, ⟨R²⟩, ΔR, ⟨1/R⟩, ⟨1/R²⟩, T_kin = E − ⟨V⟩ (including rotation), T_rad, the Lindemann ratio ΔR/⟨R⟩, and the radial density P_v(R) = |u_v(R)|² on the grid (∫P dR = 1).
5. **Critical masses** (`scan.critical_masses`): bisection in log M on the node count, run at two Numerov step sizes. The count converges as O(dx⁴), so the two results are Richardson-extrapolated. **This is the BO prediction**, without the adiabatic correction.
6. **Mass scan** (`scan.mass_scan`): writes a CSV with the following columns:
   * Thresholds (BO, exact, and their difference).
   * E0, binding and err_E0.
   * Level count (node count and DVR).
   * v=0→1 spacing.
   * E(J=1)−E(J=0) and E(J=2)−E(J=0).
   * ⟨R⟩, ΔR, Lindemann ratio, ⟨1/R²⟩, B0 = ⟨1/R²⟩/(2mu_N), T_kin.
   * R_e, E_min and k for this mu_e.
   * Harmonic estimates: ω = sqrt(k/mu_N), ΔR_harm = (2 mu_N ω)^(−1/2), E_min + ω/2.
   * Optionally, E0 with the adiabatic correction and its threshold.

## How to run (from `solvers/`)

```bash
cd solvers
PY=python3
$PY -m bo_h2plus electronic --R 1 2 5 [--state 2pu] [--mu 0.9997]
$PY -m bo_h2plus equilibrium --M 1836.15267245
$PY -m bo_h2plus levels --M 1836.15267245 --J 0 [--adiabatic] --out lv.csv --density dens.csv
$PY -m bo_h2plus scan --M 1 10 100 1836.15267245 1e4 1e5 --out scan.csv [--adiabatic]
$PY -m bo_h2plus scan --logspace 0 5 26 --out scan.csv
$PY -m bo_h2plus critical --nmax 10
$PY -m pytest -p no:cacheprovider bo_h2plus/tests -q   # 35 tests, ~3 s with cached tables
```

The same operations are available from Python: `from bo_h2plus import solve_levels, count_levels, mass_scan, critical_masses, electronic_energy`.

Timing (busy 10-core Mac): M ≤ 1e4 takes under 0.2 s per (M, J) and M = 1e5 about 5 s. M = 1e6 should work, but the dense DVR grows as M^{3/4}, so expect minutes. BLAS threads are capped at 4 (`BO_H2PLUS_THREADS`). Oversubscribed OpenBLAS was 40x slower on a loaded machine.

## Validation

"(memory)" marks a target quoted from memory by the requester; it is not an independent literature check.

| Quantity | Computed | Target | Difference |
|---|---|---|---|
| U(R=2), electron mass 1 | −0.602634214495265 | −0.6026342144949 (memory) | 4e-13 |
| R_e | 1.99719332 | ~1.9972 (memory) | – |
| E_min | −0.6026346191 | ~−0.6026346 (memory) | – |
| k = U''(R_e) | 0.1029705 | – | – |
| 2pσu U(R=2) | −0.1675343922 | −0.16753439 (memory) | < 1e-8 |
| R→0 | (E_el + 2)/R² → 8/3 | −2 (He⁺) | passes |
| R→∞ fit of (U + 1/2)R⁴ | −2.2500 and −7.5/R² | −9/4 | < 1e-4 |
| H2+ BO E(v=0, J=0) | −0.597233978550 | non-BO −0.597139063 (memory) | −9.49e-5 (BO lower, as expected) |
| H2+ BO + adiabatic E(0,0) | −0.597138541554 | −0.597139063 (memory) | +5.2e-7 |
| H2+ bound J=0 levels | 20 (v=0..19), both methods | ~20 (memory) | 0 |
| v=19 binding | 0.715 cm⁻¹ (BO); 0.735 cm⁻¹ (BO+ad) | – | – |
| H2+ v=0→1 | 2191.144 cm⁻¹ (BO); 2191.300 (BO+ad) | ~2191 (memory) | – |
| H2+ harmonic ω (mu_e) | 2323.40 cm⁻¹ | ~2322 (memory) | +1.4 (see note) |
| H2+ E(J=1)−E(J=0), E(J=2)−E(J=0) | 58.234, 174.237 cm⁻¹ (B₀ ≈ 29.12) | B ~30 (memory) | – |
| ⟨1/R²⟩/(2mu_N) for v=0 | 29.156 cm⁻¹ | – | – |
| H2+ ⟨R⟩, ΔR, ΔR/⟨R⟩ (v=0) | 2.06393, 0.23142, 0.1121 | harmonic ΔR 0.22682 | – |

The harmonic ω here is sqrt(k/mu_N) at the minimum of mu_e U(mu_e R). The "~2322" target is probably the spectroscopic ω_e, which is defined differently, so a 1.4 cm⁻¹ difference is not a discrepancy.

BO + adiabatic lies 5.2e-7 hartree (0.11 cm⁻¹) above the quoted non-BO energy. That is the expected sign and magnitude for the remaining non-adiabatic correction. I could not independently check the size of the non-adiabatic correction.

### Internal cross-checks (all in `tests/test_bo.py`)

* **Independent radial method**: direct ODE shooting (solve_ivp, DOP853) of the ξ equation agrees with the Jaffé recurrence to 1–2e-13 at R = 1, 2, 4. This also confirms the sympy-derived recurrence.
* **Truncation**: doubling both basis sizes changes E_el by < 2e-11 for R = 0.05 to 120. At R = 2 the change is 1e-12, which is round-off level.
* **Virial theorem** (Euler form 2E − Z∂E/∂Z + R∂E/∂R = 0): residual ≤ 3.5e-9, limited by the finite difference.
* **Mass scaling**: E(R; mass mu) computed as E(R; Z=mu)/mu agrees with mu·E(mu R; 1) to < 1e-12. The interpolated V_BO agrees with the direct solver to < 2e-11.
* **Interpolation table** vs direct solver: max 4e-12 (E_el) and 1.4e-9 (D, which is about 1e-12 hartree in E_ad).
* **Adiabatic correction**:
  * ⟨ψ|ψ⟩ = 1 and ⟨ψ|∂_Rψ⟩ = 3e-12.
  * Quadrature converged to 1e-13 (80 vs 160 vs 240 nodes).
  * An independent Cartesian finite-difference evaluation of ∂_Rψ at fixed ρ agrees to 1e-10 at R = 2.
  * Large-R limit: D(30) = 0.2499959 → 1/4, with radial part 1/12 and angular part 1/6 as derived analytically. So E_ad(∞) = mu_e²/(4M).
  * D(2) = 0.174788834, giving E_ad(R=2) = 9.514e-5 hartree for H2+.
* **Nuclear grid convergence**: for H2+, a finer grid (ppw 3, box ×2) changes every level by < 1e-10 and ⟨R⟩ by < 1e-6. The automatic per-level err_E is < 1e-12 for v ≤ 10 and ≤ 6e-11 for the most diffuse levels.
* **Level counting**: DVR and node count agree at M = 1, 5, 50, 1836, 1e4 and 1e5 (1, 1, 3, 20, 47, 150 levels; the M = 50 count was misprinted as 4 in an earlier version, rechecked 2026-11-13).
* **Critical masses**: Richardson error estimates are shown below. At M_c(1 ± 2e-3) the DVR finds n−1 levels below and n levels above.

### BO critical masses (J = 0; mass at which level v first becomes bound below −mu_e/2)

| n | v | M_c | m_eff = (2M+1)/4 | est. error |
|---|---|---|---|---|
| 1 | 0 | 0.964915424 | 0.732457712 | 1e-12 |
| 2 | 1 | 11.0814565936 | 5.79072830 | 2e-10 |
| 3 | 2 | 31.16777427 | 15.83388714 | 4e-9 |
| 4 | 3 | 60.8466744 | 30.6733372 | 3e-8 |
| 5 | 4 | 99.7973934 | 50.1486967 | 1e-7 |
| 6 | 5 | 147.663680 | 74.0818401 | 4e-7 |
| 7 | 6 | 204.060272 | 102.280136 | 1e-6 |
| 8 | 7 | 268.640823 | 134.570412 | 2e-6 |
| 9 | 8 | 341.228978 | 170.864489 | 5e-6 |
| 10 | 9 | 421.94099 | 211.220493 | 1e-5 |

sqrt(m_eff) at successive thresholds is spaced by about 1.55, as WKB predicts. The estimated error is |M(dx) − M(dx/2)|/15. The Richardson value should be better than that, but this has not been proven.

## Accuracy you can claim

* **E_el(R)**: ~1e-12 hartree (round-off limited) for 0.01 ≤ R ≤ 200. It has been checked by an independent method at three distances, and against the target at R=2 to 4e-13.
* **BO vibrational levels**: H2+ to ≲1e-10 hartree, with the per-level estimate attached to every result. For M ≳ 300 typically 1e-13–1e-11. For small M (≲ 30), ~1e-9 (worst seen: 1.2e-9 at M = 3). At small M the wavefunction reaches R = 0, where the Coulomb behaviour u = R(1 + mu_N R + …) limits the sinc-DVR. A finer map near the origin (R_s = 0.05) is used there.
* **Adiabatic correction**: D(R) to ~1e-9 relative, which is ~1e-12 hartree in E_ad.

## Known limitations and what was *not* verified

* There are no literature tables for the adiabatic correction or the critical masses. The adiabatic correction is validated by internal consistency, two analytic limits, and the comparison of H2+ BO+ad with the quoted (memory) non-BO energy. The critical masses are validated internally only.
* The tail coefficients −9/4 and −15/2 come from a numerical fit (−15/2 only to ±0.3). The R⁻⁷ term is not included beyond 200 bohr; its effect there is < 1e-15.
* Node counting with `adiabatic=True` uses the BO C4 plus a fitted R⁻⁴ coefficient of D. The critical masses are computed for pure BO only.
* For J>0, beyond R_b the continuation neglects the C4 tail relative to the centrifugal term. R_b is ≥ 30β, so this is negligible but not tested near a J>0 critical mass.
* Expectation values are DVR quadratures. ⟨R⟩ agrees between grids to ~1e-10–1e-7 (looser for near-threshold levels, see err_R).
* `P_v(R)` is given on the nonuniform DVR grid, and trapezoid integration of it is only good to ~1e-3. Use Σ P_i R'(x_i) dx for exact normalization.
* **BO at small M**: at M = 1 BO predicts **one** extremely weakly bound level: E = −0.33334422, bound by 1.09e-5 hartree below −mu_e/2 = −1/3, with ⟨R⟩ ≈ 154 bohr and ΔR/⟨R⟩ ≈ 0.95. This is because M = 1 lies just above the first BO critical mass, 0.9649. It is physically meaningless: the exact three-body threshold is −1/4, and the real Ps⁻-like system is bound by a completely different mechanism. With the adiabatic correction included, the M=1 "level" lies 7e-3 below its own threshold.
* 2pσu only has the electronic energy (no nuclear levels, no adiabatic correction).

## Files

* `electronic.py`: two-center solver, equilibrium, wavefunction.
* `adiabatic.py`: D(R), E_ad.
* `potential.py`: tables, mass scaling, V_BO, E_ad.
* `nuclear.py`: DVR levels, node counting, thresholds.
* `scan.py`: mass scan CSV, critical masses.
* `__main__.py`: CLI.
* `tests/test_bo.py`: 35 tests.
* `_cache/`: generated tables.

## Unequal nuclear masses (added 2026-11-16, new module `unequal.py`)

`unequal.py` (nothing else changed) covers {A+ (M), B+ (rM), e-}: mass constants and thresholds (`Masses`,
exact threshold −μ_Be/2), the exact linear map λ = M_eff/μ_e = 2rM/(1+r) + 2r/(1+r)² (`lam_of_M`, `M_of_lam`),
the electron kinetic ⟨p²⟩ from the virial theorem, and the 1sσg/2pσu coupling functions (D_u, t_g, t_u,
⟨g|∂_z|u⟩, k_gu; `gu_functions`, cached table `_cache/gu_v1.npz`, rebuilt in ~30 s). Theory and use:
`paper_data/unequal_masses/THEORY.md`. Tests: `tests/test_unequal.py` (9 tests: hypervirial, virial, asymptotic
atomic reduced-mass eigenvalues, DBOC vs Schwenke 2001 HD+ table to 1e-3).
