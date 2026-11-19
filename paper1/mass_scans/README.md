# Exact versus Born–Oppenheimer mass scans of {X⁺, X⁺, e⁻}

*Last updated: 2026-11-17 20:42*

> **Mass constants (2026-11-17):** Paper 1 quotes physical-mass values at CODATA 2022 (m_p/m_e = 1836.152673426, m_μ/m_e = 206.7682827); they are in `paper1/codata2022/` (README, results.csv) and `solvers/threebody/results/codata2022/`. Examples and tables in this README that name other masses (CODATA 2018 1836.15267343, 206.768283; CODATA 2010 1836.15267245, 206.7682843, used for the King et al. comparison and called "CODATA 2012" in older notes) are kept as documentation of those runs; they differ from the CODATA 2022 values only in the last digits.

Written 2026-11-13. Two identical particles of mass M (charge +1) and one of mass 1 (charge −1), atomic units, κ = M^(−1/4). Exact (non-BO) results come from `solvers/threebody`; BO and BO + adiabatic (BO+ad) results from `solvers/bo_h2plus`. Both packages were used unchanged; everything here is new scripts and output. This answers referee concerns M9 (an exact BO baseline instead of harmonic estimates) and M11 (internal-coordinate descriptors first, the CM density second and labelled as descriptor-dependent).

## How to regenerate

```bash
cd paper1/mass_scans
export THREEBODY_BLAS_THREADS=3 BO_H2PLUS_THREADS=3
PY=python3
$PY scan.py bo        # BO and BO+ad, J = 0, 1, 2, derivatives, P(R)  (~20 s)
$PY scan.py exact     # exact L = 0 'S', roots 0 and 1, all descriptors (~35 min on a loaded 10-core Mac)
$PY scan.py rot       # exact L = 1 'A' and L = 2 'S' lowest levels (~45 min)
$PY scan.py exact 0.62 0.64 0.66 0.68 1.05 1.1 1.15 ; $PY scan.py bo 0.62 0.64 0.66 0.68 1.05 1.1 1.15
$PY analyze.py        # fits.txt, energies/descriptors/spacings.csv, fig1..fig6
```
Every stage appends one row per mass and skips masses already present, so an interrupted run is resumed by repeating the command. `common.py` holds the mass grid and the CM-density routine. Wall time used here: about 65 min with two or three processes of 3 BLAS threads each.

## Files

| file | content |
|---|---|
| `scan.py`, `analyze.py`, `common.py` | compute, analyse, shared settings |
| `exact_L0.csv` | per M: basis, scales, energy sequences (`E0_seq`, `E1_seq`), extrapolations, all descriptors at the largest basis and (`*_prevN`) at the next-smaller one |
| `exact_rot.csv` | per M and L: lowest L = 1 'A' and L = 2 'S' level, sequences, extrapolation, ⟨R⟩ |
| `bo.csv` | BO and BO+ad: E(v,J) for v = 0, 1 and J = 0, 1, 2, ⟨R⟩, ΔR, Lindemann, ⟨1/R²⟩, ⟨T_R⟩, dE/dlnM, d²E/d(lnM)² |
| `energies.csv` | exact E0 (variational and extrapolated, error), E_BO, E_BO+ad, the differences, E/μ_Xe |
| `descriptors.csv` | exact descriptors with basis-change errors, and the BO / BO+ad counterparts |
| `spacings.csv` | vibrational (v=1 − v=0) and rotational (L=1, L=2 − L=0) spacings, exact / BO / BO+ad, ratios, relative errors |
| `gR_exact.csv`, `gR_bo.csv` | X–X distance distribution P(R) (= 4πR² g_XX) on R ∈ [0, 20] bohr and on R/⟨R⟩ ∈ [0, 4] |
| `density_cm.csv` | one-particle density of X about the centre of mass, ρ_X(s), s/ℓ ∈ [0, 3] |
| `fits.txt`, `density_crossover` block in it | all fits, local slopes, crossover estimates |
| `fig1_energies_bo_error.png` | E vs M; \|E − E_BO\| and \|E − E_BO+ad\| vs κ (log–log) with κ⁴, κ⁶ guides; E/μ_Xe (H6) |
| `fig2_lindemann_R.png` | Lindemann ratio and ⟨R⟩ vs M, exact / BO / BO+ad |
| `fig3_gR.png` | P(R) vs R/⟨R⟩ for selected M; unscaled; exact − BO difference (M = 20, 206.8) |
| `fig4_density_cm.png` | ρ_X(s)/ρ_X(0) about the CM and its curvature at s = 0 vs M, with literature crossovers |
| `fig5_mass_derivatives.png` | −dE/dlnM and d²E/d(lnM)², exact / BO / BO+ad |
| `fig6_spacings.png` | spacings, d₂/d₁ and vib/d₁, relative errors of BO and BO+ad |
| `data/L0_M*.npz` | exact wave functions (roots 0, 1) at every M (threebody `Result.save` format) |
| `run_log.txt`, `log_*.txt`, `PROGRESS.md` | logs and progress notes |

## Method

**Mass grid.** 30 masses from 0.05 to 10⁴ (0.05, 0.1, 0.2, …, 1, 1.2, 1.5, 2, 3, 4.5, 6, 9.5256, 12, 15, 20, 29.639, 50, 100, 206.768283, 500, 1000, 1836.15267343, 4000, 10⁴), plus 0.62–0.68 and 1.05–1.15 added to locate two crossovers. **Mass constants:** the scan grid (and every CSV, npz and figure in this folder) used CODATA 2018 at the two physical points (206.768283, 1836.15267343); the physical-mass values quoted in Paper 1 and in the two physical-mass table rows of this README are CODATA 2012 (m_μ = 206.7682843, m_p = 1836.15267245, as King et al. 2013), recomputed in `paper1/codata2012/` (`results.csv`, `changes.md`). The difference is below the printed precision except in the last digits listed in `changes.md` (Mu₂⁺ E₀ −6.0e−11, H₂⁺ E₀ +1.6e−12). The CSV files here are unchanged and still hold the CODATA-2018 rows. 4.5 and 12 sit just above the L = 1 and L = 2 critical masses (4.208, 11.489), 9.5256 and 29.639 at the v = 1 and v = 2 critical masses.

**Exact L = 0.** Perimetric Laguerre basis, exchange-symmetric. Scales (α_X, α_e) optimized for root 0 by Nelder–Mead at N_opt (12, total degree, for M < 150; 20 with weights (1,1,2) above), warm-started from the neighbouring mass, then **fixed**; nested sequence N = 16, 20, 24, 28 (M < 150, size 2360), 28–40 (150 ≤ M < 3000, size 3311), 32–44 (M = 4000, size 4324), 36–48 (M = 10⁴, size 5525). Energies are extrapolated with the package's geometric-tail rule (heuristic error). Root 1 (v = 1) is taken from the same basis and, for 9.5 ≤ M < 150, also from a separate basis with scales optimized for root 1 (the lower of the two upper bounds is used).

**Descriptors** (root 0, largest basis): ⟨R⟩, ⟨R²⟩, Lindemann ratio ΔR/⟨R⟩, ⟨1/R⟩, ⟨1/R²⟩ (R = X–X distance), T, V, virial ratio, the Hellmann–Feynman X kinetic energy K_X = (⟨T₀⟩ + ⟨T₁⟩)/M = −dE/dlnM (exact for the variational eigenvalue at fixed basis), and d²E/d(lnM)² as the central difference of the Hellmann–Feynman derivative at M e^(±h), same matrices, h = 0.02 and 0.01, Richardson-combined. Its error = |Richardson correction| + |change from the N−4 basis|. Checks stored in the CSV: the HF first derivative agrees with a finite difference of E to O(h²), and the second derivative from the energy alone agrees to ~1e−6 relative or better. Every descriptor is also evaluated at the second-largest basis (`*_prevN`) as convergence evidence.

**BO and BO+ad.** `bo_h2plus.solve_levels` for J = 0, 1, 2 with and without the diagonal adiabatic correction; each model's levels are counted at its own threshold. dE/dlnM and d²E/d(lnM)² by central differences (h = 0.02, 0.01, Richardson). BO has no bound level below M = 0.9649 (its own critical mass), so BO columns are empty there; BO+ad is bound at all M tested.

**Rotational levels.** `threebody.rotational` (natural parity; L = 1 antisymmetric and L = 2 symmetric are the rotational tower of the ground manifold). Scales optimized for L = 1 root 0 (N_opt 12 or 16 with (1,1,2); 20 at M = 10⁴), reused for L = 2; sequences N = 14, 18, 21 (M < 150), 24, 28, 32 (M < 3000), 26, 30, 34 (M = 4000; L = 2 at 10⁴), 30, 34, 38 (L = 1 at 10⁴). Computed only above the exact critical masses.

**Density of X about the centre of mass** (not provided by the solver; computed here, `common.density_about_cm`). Mathematically exact reduction: with R = x₀ − x₁ and r = x₂ − (x₀ + x₁)/2, the vector from the CM to X₀ is s = R/2 − r/(2M+1). Changing variables (R, r) → (s, r) (Jacobian 8) gives ρ_X(s) = 8 ∫ d³r |Ψ(R = 2s ẑ + 2r/(2M+1), r)|², independent of the direction ẑ because L = 0. Ψ depends only on the three distances, which are elementary functions of (s, r, θ), so the remaining 2D integral is done by Gauss–Legendre quadrature (log-spaced panels in r, θ rather than cos θ so that the coalescence kinks sit at endpoints). Checks: ∫4πs²ρ ds = 1 to 3e−5 and ⟨s²⟩ against the exact value from ⟨r₀₁²⟩, ⟨r₀₂²⟩ to 3e−4 relative at M = 0.3, 1, 20 (both limited by the outer s-quadrature of the check, not by ρ); doubling the inner quadrature changes ρ by < 2e−6. The topology indicator is the curvature a₂ in ρ(s)/ρ(0) = 1 + a₂(s/ℓ)² + a₄(s/ℓ)⁴ fitted on s/ℓ ≤ w, ℓ = ⟨s²⟩^(1/2), w = 0.05, 0.08, 0.1.

## Convergence evidence

* **Ground-state energies** (`energies.csv`, column `E0_seq`): every sequence decreases monotonically to the rounding level except where flagged. Last increments: M = 0.05: 1.5e−10, 3.1e−11, 8.4e−12; M = 1: 9.6e−10, 2.1e−10, 5.8e−11 (Ps⁻ −0.262005070224 ± 2e−11, as in the package validation); M = 9.5–1000: < 1e−12 (rounding floor); M = 1836.15: 6.9e−11, 1.4e−12, 1.2e−13 (scan-grid E = −0.597139063081 at CODATA 2018; at the CODATA-2012 mass E = −0.5971390630793, within 1e−13 of the verified reference −0.597139063079 for m_p = 1836.15267247); M = 4000: 7.6e−10, 4.8e−11, 1.7e−12; **M = 10⁴: 6.3e−8, 1.5e−8, 2.0e−9, extrapolated error 2.2e−9** (the only mass where the energy error is not negligible: 5% of E − E_BO+ad there). Estimated errors: ≤ 3e−11 for M ≤ 1, ≤ 1e−12 for 3 ≤ M ≤ 4000.
* **Descriptors** (change from the N−4 basis, `descriptors.csv` `*_basis_err`): ⟨R⟩ ≤ 3e−9 relative for M ≤ 1 (2.4e−7 bohr absolute at M = 0.05, where ⟨R⟩ = 94 bohr), < 5e−10 bohr for 3 ≤ M ≤ 4000, 3.7e−7 bohr at 10⁴. Lindemann ratio ≤ 3e−9 (M ≤ 1), < 5e−10 for 3 ≤ M ≤ 4000, 5e−7 at 10⁴. dE/dlnM ≤ 1.4e−8 (at 10⁴; ≤ 4e−10 elsewhere). d²E/d(lnM)²: finite-difference part ≤ 9e−7 hartree, basis part ≤ 2e−9 for M ≤ 4000 and 9e−8 at 10⁴. Virial ratio −V/T = 2 to ≤ 4e−9.
* **v = 1**: converged to ≤ 2e−10 for M ≥ 12 (root-1 scales) and 1e−13…5e−11 for 50 ≤ M ≤ 4000; 3e−8 at 10⁴. At M = 9.5256 the state sits at its threshold and is not bound in any basis (as it must be).
* **Rotational levels**: extrapolated errors ≤ 1e−12 for 6 ≤ M ≤ 1836 (L = 1 at 4.5: 7e−10; L = 2 at 12: 5e−8, close to its critical mass), 1.5e−10 at M = 4000, and at 10⁴, after the 2026-11-14 extension (`paper1/convergence_large_M/`), 8e−13 (L = 0), 1.4e−10 (L = 1) and 1.4e−8 (L = 2), and 5.0e−9 for L = 2 after the N = 46 run of 2026-11-15; d₂/d₁ = 2.99853 ± 0.00010, below 3 and equal to BO (2.99853) like every other mass. The earlier ±0.002 value 3.0009 was an extrapolation artefact of the short N ≤ 34 sequence. fits.txt was refitted 2026-11-14 with these values (`analyze.py` reads `paper1/convergence_large_M/results.csv` and overrides the 10⁴ row). The changes are small: E − E_BO+ad exponent 6.070 ± 0.013 → 6.067 ± 0.014, two-term a = −0.039339, d₁ exponent unchanged at 3.952 (window M ≥ 1000). The H₂⁺ (0,1) level −0.5968737388 and (0,2) level −0.5963452055 (CODATA-2012 mass; the scan-grid CODATA-2018 values were −0.5968737387851 and −0.5963452054905; the CODATA-2012 digits beyond the 10th are at the extrapolation error 4.5e−11 and are not quoted) differ from the verified references (paper1/rotational, m_p = 1836.15267247) by −4e−13 and −1.4e−13; the mass difference alone shifts them by ≈ −1.5e−12, so the agreement is at the 1e−12 level. At 10⁴ the first L = 1 run (N ≤ 34) had an error of 5e−7; it was rerun with N_opt = 20 and N = 30, 34, 38 (error 2.1e−8, d₁ known to 5e−4 relative). The L = 2 level at 10⁴ (N ≤ 34, error 1e−7) is not converged enough for d₂/d₁ (3.0009 ± 0.002, above the rigid-rotor 3, unlike every other mass).
* **BO side**: per-level DVR error ≤ 1e−9 (small M) to 1e−12; adiabatic correction to ~1e−12 hartree (package README).

## Results

### 1. Energies and the BO error (fig. 1, `energies.csv`, `fits.txt`)

*The rows at 206.768 and 1836.15 below show the CODATA-2012 physical-mass values (see the Mass-constants note under Mass grid); `energies.csv` holds the CODATA-2018 scan values.*

| M | E exact (extrap.) | err | E_BO | E_BO+ad | E − E_BO | E − E_BO+ad |
|---|---|---|---|---|---|---|
| 0.05 | −0.025062727903 | 3e−12 | (unbound) | −0.0041766628 | – | −2.09e−2 |
| 1 | −0.262005070224 | 2e−11 | −0.3333442151 | −0.2293938855 | +7.13e−2 | −3.26e−2 |
| 3 | −0.400267038840 | 4e−12 | −0.4408310610 | −0.3899510024 | +4.06e−2 | −1.03e−2 |
| 20 | −0.536021015599 | 2e−14 | −0.5441608059 | −0.5354224604 | +8.14e−3 | −5.99e−4 |
| 206.768 | −0.585126098252 | 3e−13 | −0.5859615137 | −0.5851111765 | +8.354e−4 | −1.492e−5 |
| 1000 | −0.595093299696 | 4e−13 | −0.5952673546 | −0.5950919829 | +1.7405e−4 | −1.317e−6 |
| 1836.15 | −0.597139063079 | 1e−13 | −0.5972339785 | −0.5971385416 | +9.4915e−5 | −5.215e−7 |
| 4000 | −0.598954614841 | 1e−12 | −0.5989982307 | −0.5989544548 | +4.3616e−5 | −1.600e−7 |
| 10⁴ | −0.600328524648 | 8e−13 | −0.6003459839 | −0.6003284846 | +1.7459e−5 | −4.005e−8 (±0.005e−8) |

Signs: BO lies below the exact energy at every M where it is bound; BO+ad lies above it at every M (both expected: the adiabatic correction is positive, the remaining non-adiabatic correction is a second-order lowering).

**Fitted exponents (free fits of ln|ΔE| vs ln κ; window M ≥ 1000, error = max(statistical, spread over windows M ≥ 500, 1000, 1836)):**
* **E − E_BO ∝ κ^p with p = 3.995 ± 0.002**; local slopes rise monotonically 3.80 (M 20–30) → 3.997 (M 4000–10⁴). Prefactor: M(E − E_BO) = 0.1628 (M = 20), 0.1737 (500), 0.17428 (1836), 0.17459 (10⁴); a fixed-exponent fit a κ⁴ + b κ⁶ on M ≥ 500 gives a = 0.17498, b = −0.030. Theory: E − E_BO = ⟨E_ad⟩ + O(κ⁶) with ⟨E_ad⟩ → D(R_e)/M, and D(2.0) = 0.1748 in bo_h2plus. **p = 4 confirmed, coefficient consistent with the adiabatic correction at R_e.**
* **E − E_BO+ad ∝ κ^p with p = 6.07 ± 0.01 (free fit), local slopes decreasing monotonically 6.38 → 6.056 (last pair; with the 10⁴ error, ±0.2 on that slope)**, i.e. approaching 6 from above with a κ⁸ correction of the same sign. M^(3/2)(E − E_BO+ad) = −0.0535 (M = 20), −0.0426 (500), −0.0410 (1836), −0.0400 (10⁴); fixed-exponent fit a κ⁶ + b κ⁸ (M ≥ 500): a = −0.0393, b = −0.073. Theory (justified, leading order): the first non-adiabatic energy correction is second order in the coupling (1/M)⟨φ_n|∂_Rφ₀⟩∂_R; with ⟨(∂_Rχ)²⟩ ∝ μ_N ω ∝ M^(1/2) this gives O(M^(−3/2)) = κ⁶, negative; equivalently, a correction of relative order 1/M to the effective vibrational mass acting on the zero-point energy ω/2 ∝ κ². **So the κ⁶ law is consistent with the data, the free exponent is 6.07 ± 0.01 because the window is not fully asymptotic; I would quote "κ⁶ with coefficient −0.040 ± 0.001 (asymptotic)".** The coefficient's uncertainty is my judgement from the drift of M^(3/2)ΔE, not a fit error.
* **Normalized energy E/μ_Xe (H6):** non-monotonic, with a maximum −0.52256 at M ≈ 0.5 on the grid (−0.5263 at M = 0.05, −0.6004 at 10⁴). H6 is confirmed in the sense stated (a maximum near M ~ 1); the maximum is shallow and its position is only grid-resolved (0.4–0.6).

### 2. Descriptors (figs. 2, 3, 5; `descriptors.csv`)

| M | ⟨R⟩ | ΔR/⟨R⟩ exact | BO | BO+ad | −dE/dlnM | d²E/d(lnM)² |
|---|---|---|---|---|---|---|
| 0.05 | 93.667 | 0.54584 | – | 0.908 | 0.02381 | −0.02144 |
| 0.3 | 19.679 | 0.55066 | – | 0.785 | 0.09247 | −0.04988 |
| 1 | 8.5486 | 0.52446 | 0.970 | 0.620 | 0.13324 | −0.00608 |
| 3 | 4.7904 | 0.45323 | 0.522 | 0.476 | 0.10971 | +0.04061 |
| 9.5256 | 3.2407 | 0.36663 | 0.37661 | 0.37057 | 0.06203 | +0.03590 |
| 20 | 2.7754 | 0.31605 | 0.31924 | 0.31743 | 0.03968 | +0.02451 |
| 206.768 | 2.20522 | 0.18895 | 0.18908 | 0.18902 | 0.009579 | +0.005573 |
| 1836.15 | 2.06391 | 0.11212 | 0.11212 | 0.11212 | 0.002847 | +0.001521 |
| 10⁴ | 2.02541 | 0.07392 | 0.07392 | 0.07392 | 0.001171 | +0.000604 |

* **Lindemann ratio.** Large-M exponent (exact) **x = 0.980 ± 0.004** in κ (window M ≥ 1000); local slopes 0.83 (M 20–30) → 0.986 (4000–10⁴), so the asymptotic value is 1 (BO harmonic ΔR ∝ κ), approached slowly from below. BO gives the same numbers (0.980) at large M. Prefactor ≈ 0.71 (window-dependent, drifting up). At small M the exact ratio is almost constant, 0.52–0.55, with a shallow **maximum 0.5507 at M ≈ 0.2–0.3**, below the two-body 1s value 1/√3 = 0.577; BO and BO+ad overestimate it strongly for M ≲ 3 (BO+ad 0.62 at M = 1, 0.91 at 0.05). The exact − BO Lindemann difference scales with local slopes ≈ 5.02–5.11 for 207 ≤ M ≤ 4000 (consistent with κ⁵ = κ × relative κ⁴; the last point is affected by the 10⁴ basis error).
* **⟨R⟩ − R_e** (R_e = 1.99719, BO): exponent 2.036 ± 0.008 (local slopes → 2.025), i.e. κ² (anharmonic shift ∝ ΔR²). ⟨R⟩_exact − ⟨R⟩_BO: local slopes 6.5 → 5.76, not yet asymptotic; I do not quote an exponent.
* **X kinetic energy −dE/dlnM**: exponent 2.11 ± 0.02 (local slopes → 2.08), consistent with ω/4 ∝ κ² plus a κ⁴ correction. Non-monotonic in M: maximum 0.13355 at **M ≈ 1.1**; d²E/d(lnM)² changes sign there (exact: 1.110, linear interpolation between 1.1 and 1.15; BO+ad: ≈ 1.10). Large-M exponent of d²E/d(lnM)²: 2.20 ± 0.04 (local slopes → 2.15), consistent with ω/8 ∝ κ² plus corrections. BO+ad reproduces both derivatives to ≤ 1% for M ≥ 50 and deviates by tens of percent for M ≲ 3. The exact first derivative has no feature other than this smooth maximum (H1: no sign of a non-analyticity on this grid).
* **Pair distribution** (fig. 3): P(R) vs R/⟨R⟩ is broad and almost M-independent for M ≤ 1 (the curves for 0.1, 0.5, 1 nearly coincide) and sharpens continuously above M ~ 3. The exact − BO difference of P(R) (panel c) is a smooth shift towards larger R for M = 20 and 206.8; at M ≳ 1000 the BO P(R) interpolated from its DVR grid is too noisy to show the difference (not shown).
* **Coalescence** ⟨R⟩³ g_XX(0): 0.228 (M = 0.05) → 0.107 (1) → 1.5e−8 (50). ln g(0) vs √M = κ^(−2): local slopes −1.65 (M ≈ 1) → −2.93 (M = 30–50), still drifting, so the expected ln g(0) ∝ −c κ^(−2) is **not** established (the slope is not constant); beyond M ≈ 50 g(0) falls below the basis resolution (values non-monotonic, noise).

### 3. Density of X about the centre of mass (fig. 4; descriptor-dependent, reference point = CM)

* The curvature a₂ changes sign once, from a maximum at the CM (a₂ < 0) to a hole (a₂ > 0), between M = 0.62 (a₂ = −0.45) and 0.64 (a₂ = +0.12): **M_c(CM density) = 0.635 ± 0.01**, where the error is the spread over fit windows w = 0.05–0.1 (0.627–0.644); the basis change of a₂ (N vs N−4) is ≤ 7e−3 near M_c, which moves M_c by < 1e−3. a₂ is smooth and monotonic across M_c (no sign of a singularity). Above M_c the maximum moves out of the origin; ρ(0)/ρ_max = 0.92 at M = 1, 0.32 at M = 3, 4e−4 at M = 20.
* **Comparison with the literature (RESULTS_SUMMARY §2):** our 0.635 lies inside the range 0.4–0.8 attributed to Mátyus et al. 2011 (from secondary sources; not checked against the paper). King et al. 2013 (0.491) and Salas et al. 2022 (≈ 0.98) use **different** criteria (binding-ratio minimum, marginal-distribution curvature), so they are not expected to coincide; the spread 0.49–0.98 between criteria is larger than any of our numerical uncertainties, which is the descriptor dependence the referee (M11) and Ludeña et al. warned about.
* Other "crossover" markers from the same wave functions, for H5: Lindemann maximum at M ≈ 0.2–0.3; E/μ_Xe maximum at M ≈ 0.5; CM-density topology at 0.635; maximum of the X kinetic energy (zero of d²E/d(lnM)²) at 1.11; first rotational level bound at 4.208, first vibrational excitation at 9.526. These are distinct, and none is sharp.
* Limitations: at small M the CM approaches the electron and ρ develops a near-cusp, so |a₂| becomes large (−59 at M = 0.05); only the sign is used. For M ≳ 50, ρ(0) is < 1e−7 of the maximum and its curvature is basis noise (a₂ even turns negative at M ≥ 100 with N-to-N−4 changes of order 100%); those values are not physics and are excluded from fig. 4b. About the pair midpoint the density is 8 g_XX(2s), which has a minimum at s = 0 at every M (positive Kato cusp of the repulsive X–X pair), so that reference point gives no transition at all.

### 4. Vibrational and rotational spacings (fig. 6, `spacings.csv`)

| M | E(v=1)−E(v=0) | BO+ad rel. err. | d₁ = E(L=1)−E(L=0) | M·d₁ | 2⟨1/R²⟩ | d₂/d₁ exact | BO | vib/d₁ |
|---|---|---|---|---|---|---|---|---|
| 4.5 | (unbound) | – | 3.1109e−2 | 0.1400 | 0.227 | (L=2 unbound) | – | – |
| 12 | 4.9710e−2 | −5.4e−3 | 2.2262e−2 | 0.2671 | 0.325 | 2.218 | (unbound) | 2.23 |
| 20 | 5.1166e−2 | +3.2e−3 | 1.6120e−2 | 0.3224 | 0.364 | 2.4747 | 2.4638 | 3.17 |
| 50 | 4.3212e−2 | +2.7e−3 | 7.8985e−3 | 0.3949 | 0.415 | 2.7489 | 2.7476 | 5.47 |
| 206.768 | 2.62977e−2 | +7.6e−4 | 2.19765e−3 | 0.4544 | 0.460 | 2.9316 | 2.9315 | 11.97 |
| 1000 | 1.324075e−2 | +1.7e−4 | 4.81722e−4 | 0.4817 | 0.4829 | 2.9854 | 2.9854 | 27.5 |
| 1836.15 | 9.98338e−3 | +9.2e−5 | 2.65324e−4 | 0.4872 | 0.4878 | 2.9920 | 2.9920 | 37.6 |
| 4000 | 6.89551e−3 | +4.2e−5 | 1.22988e−4 | 0.4920 | 0.4923 | 2.9963 | 2.9963 | 56.1 |
| 10⁴ | 4.42579e−3 (±3e−8) | (unresolved) | 4.95520e−5 (±1.4e−10) | 0.4955 | 0.4956 | 2.99853 (±0.00010) | 2.99853 | 89.3 |

(H₂⁺: vib 9.98338e−3 Eh = 2191.10 cm⁻¹, d₁ = 58.232 cm⁻¹.)
* Large-M exponents of the exact spacings: vibrational **1.905 ± 0.023** (local slopes 0.50 at M = 20–30 → 1.94 at 4000–10⁴; expected 2, approached slowly from below); rotational d₁ **3.952 ± 0.013** (local slopes → 3.97 at 4000–10⁴; expected 4). Both are consistent with the BO orders κ² and κ⁴; the free fits are biased low by the next order.
* The relative error of BO+ad spacings is ∝ κ⁴ (≈ 1/M): vibrational +0.15/M (exponent 3.94 ± 0.01), rotational d₁ ≈ +0.07–0.08/M (exponent 4.04 ± 0.07). Plain BO: d₁ error changes sign between M = 100 and 207, vibrational error between 500 and 1000 (the κ⁴ adiabatic and the mass-scaling errors cancel there), so no power law is fitted.
* d₂/d₁ approaches the rigid-rotor 3 from below and BO reproduces it to 5e−4 relative at M = 50, ≤ 1e−4 for M ≥ 100 (0.4% at M = 20): the deviation from 3 is centrifugal distortion, already present in BO. M·d₁ approaches 2⟨1/R²⟩ (rigid rotor with the exact ⟨1/R²⟩): 38% below it at M = 4.5, 1% at M = 207, 0.1% at M = 1836.
* vib/d₁ grows ∝ κ^(−2) at large M (H8 in its inverted form: d₂/vib ∝ κ²); below M ≈ 12 there is no tower at all (L = 2 unbound) and at M = 12 the vibrational spacing is only 2.2 d₁: rotation and vibration are not separated.

## What is rigorous, what is fitted

* **Rigorous / converged** (to the stated errors): all exact energies at M ≤ 4000 (≤ 3e−11), the BO and BO+ad energies, E − E_BO and E − E_BO+ad at every M ≤ 4000 (differences ≥ 1.6e−7, errors ≤ 1e−12 for M ≥ 3), all descriptors except the CM curvature at M ≳ 50 and g(0) at M ≳ 50, spacings for 6 ≤ M ≤ 4000. Hellmann–Feynman derivatives are exact for the fixed-basis eigenvalue; their basis error is estimated, not bounded.
* **Less converged:** M = 10⁴ energy (±2.2e−9, i.e. 5% of E − E_BO+ad), its v = 1 and rotational levels (v = 1: ±3e−8; L = 1: ±2e−8; L = 2: ±1e−7), and therefore the last local slope of every BO-error fit.
* **Fitted / extrapolated:** all exponents (free power laws on a finite window; the quoted uncertainties are statistical plus window spread and do not include the systematic bias from the next order, which is visible in the drifting local slopes), the asymptotic coefficients (0.1750 for κ⁴, −0.040 for κ⁶), M_c = 0.635 of the CM-density topology (linear interpolation in ln M of a fitted curvature), the zero of d²E/d(lnM)² at 1.11.
* **Theory statements** (κ⁴ for BO, κ⁶ for BO+ad, κ, κ², κ⁴ for ΔR, ω, B): standard order counting in the BO expansion, sketched above; I did not derive the κ⁶ coefficient.

## Not done / open

* The 10⁴ point would need N ≥ 52 (weights 1,1,2) or a better-localized X–X basis to reach 1e−10.
* No literature comparison of E − E_BO+ad for H₂⁺ beyond the package README (5.2e−7, matching our −5.215e−7).
* g_XX(0) asymptotics and the CM density at large M need a basis that represents exponentially small amplitudes (or a WKB argument).
* The second root used for v = 1 near threshold (M = 9.5–12) uses compact bases only; the dedicated threshold bases of `paper1/critical_masses` were not used here.
