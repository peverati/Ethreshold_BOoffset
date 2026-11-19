# Fidelity susceptibilities of the exact {X⁺, X⁺, e⁻} ground state versus M

*Last updated: 2026-11-13 19:36*

Written 2026-11-13. Two particles of mass M (charge +1), one of mass 1 (charge −1), total L = 0 ground state, atomic units. Question: does any of the smooth crossovers found in `paper_data/mass_scans` (markers at M ≈ 0.2–0.3, 0.4924, 0.5, 0.543, 0.635, 1.11) hide a sharp feature that a fidelity susceptibility would reveal? Solvers used unchanged (`solvers/threebody`); all new code is in this folder.

## Results in brief

* **No susceptibility has a sharp feature anywhere on 0.05 ≤ M ≤ 10⁴, and none has a peak at any of the marker masses (0.2–0.3, 0.4924, 0.5, 0.543, 0.635, 1.11).** On the dense windows (fine ln M grids, spacing 0.02–0.08, 0.17 ≤ M ≤ 48.5) every curve is fitted by a degree-6 polynomial in ln M to residuals of 1e−9 to 2e−7, which is at or below the N/N−4 basis change. Windows built on different fixed bases agree with each other and with the per-mass-optimized grid to ≤ 3e−5 relative (typically 1e−7 to 1e−6).
* **Rigorous statement (not numerical):** the ground state is an isolated, nondegenerate eigenvalue below the threshold (gap > 0 at every M), and H(M) = (T₀+T₁)/M + T₂ + V is an analytic family in 1/M. Kato's analytic perturbation theory then makes the eigenprojection, and hence ψ(M) in fixed coordinates, real-analytic in M. χ_F is analytic; χ_inv and χ_sep are ratios of analytic functions with non-vanishing denominators (‖Dψ‖² > 0, Gram determinant > 0), so they are analytic too. Sharp features are therefore excluded for the exact problem; the numerics only show *where* the smooth ones are. (Excited-state thresholds such as v = 1 at 9.5256 do not break this: the ground state stays isolated.)
* **χ_F** (fixed coordinates) decreases monotonically from 1.220 (M = 0.05) to 0.0331 (10⁴). At small M it is a pure dilation: χ_F(M+1)²/‖Dψ_∞‖² = 0.984, 0.976 at M = 0.05, 0.1 → 1, with ‖Dψ_∞‖² = 1.36755 (infinite-mass H⁻-like state). At M = 1, 97.9% of χ_F is the dilation component (χ_F = 0.4698, χ_inv = 0.00998). This is the referee's point (M4) made quantitative: χ_F mostly measures the change of length unit.
* **χ_inv** (global dilation removed) rises from 2.4e−4 (M = 0.05) to 0.0306 (10⁴): steepest rise (inflection in ln M) at **M = 1.3856**, then a broad, shallow **maximum 0.026527 at M = 6.8945** and a shallow **minimum 0.025850 at M = 25.066**, then a slow approach to 1/32 from below. Peak prominence above the minimum: 6.8e−4 (2.6%); half-width ≈ 0.4 in ln M (factor ≈ 1.5 in M); curvature −3.8e−3 per (ln M)². Location error from N vs N−4 and spline vs local quintic: < 1e−4 in M.
* **χ_sep** (R and ρ dilations removed) is monotone increasing everywhere, from 1.24e−4 to 0.0306. Smooth inflections (in ln M) at **M = 0.2442** (a shoulder where the slope nearly vanishes, 3.1e−4 per ln M; this is inside the Lindemann-maximum range 0.2–0.3), **1.428** (steepest rise), **9.78** and **42.9** (weak). None of the other markers.
* **Large M.** Derived expectation: in the BO regime ψ ≈ φ_el(ρ; R) u(R)/R with u a Gaussian in R − R_e of width ∝ M^(−1/4), and only the width changes at leading order. For a 1D Gaussian, F = sech^(1/2)(Δ ln σ) gives χ = ½ (d ln σ/d ln M)² = 1/32. Neither dilation can absorb a width change at fixed R_e (removing it costs a shift of R_e of order R_e/σ ∝ M^(1/4)), the shift of ⟨R⟩ (∝ κ²) and electronic changes are subleading, so **all three χ → 1/32 with O(κ²) corrections**. Measured (fit c₀ + c₁κ² + c₂κ⁴, M ≥ 500): c₀ = 0.031246 (χ_F), 0.031242 (χ_inv), 0.031246 (χ_sep), against 1/32 = 0.031250; c₁ = +0.182, −0.067, −0.067. Local log-slopes of |χ − 1/32| in κ reach 2.07 (χ_F), 1.88 (χ_inv), 1.94 (χ_sep) at 4000–10⁴, approaching 2. So Coulomb behaves like harmonic model B (1/32, M-independent bond length) and not like model A (χ_inv → 3/64, χ_sep → 0). Unlike model B with g = 1, χ_inv is not monotone (the broad maximum near 6.9).
* **Small M:** χ_inv, χ_sep → 0, consistent with a pure dilation. Local log-slopes vs M between 0.05 and 0.1 are 1.53 and 1.41, still increasing as M decreases. The expected M² (mass polarization ∝ M) is plausible but **not established** on this grid.

| M | χ_F | χ_inv | χ_sep | |N − (N−4)| (χ_inv) |
|---|---|---|---|---|
| 0.05 | 1.220183 | 2.3792e−4 | 1.2442e−4 | 1e−8 |
| 0.3 | 0.809472 | 2.2673e−3 | 7.4174e−4 | 9e−9 |
| 0.5 | 0.654405 | 3.8471e−3 | 1.3195e−3 | 7e−9 |
| 0.64 | 0.583854 | 5.3351e−3 | 2.2106e−3 | 5e−9 |
| 1 | 0.469826 | 9.9814e−3 | 5.6107e−3 | 3e−9 |
| 3 | 0.263018 | 2.4006e−2 | 1.7067e−2 | <1e−12 |
| 6 | 0.177737 | 2.6485e−2 | 1.9892e−2 | <1e−12 |
| 20 | 0.095585 | 2.5886e−2 | 2.2146e−2 | <1e−12 |
| 206.77 | 0.046247 | 2.7751e−2 | 2.7183e−2 | 1e−13 |
| 1836.15 | 0.035763 | 2.9820e−2 | 2.9747e−2 | 8e−11 |
| 10⁴ | 0.033118 | 3.0601e−2 | 3.0587e−2 | **8e−6** |

Full table: `fits.txt`, `results.csv`; features: `features.txt`.

**Convergence.** N vs N−4 changes: χ_F ≤ 3e−8, χ_inv and χ_sep ≤ 1e−8 for M ≤ 1.2 (relative ≤ 4e−5 at M = 0.05, where χ_inv is small); ≤ 1e−12 for 3 ≤ M ≤ 1000; 8e−11 at 1836; 1.2e−8 at 4000; **8e−6 at 10⁴** (2.6e−4 relative; the least converged mass, as for the energies). FD-vs-perturbative χ_F: ≤ 4e−9 relative at every mass. Uniform-scaling test (χ_inv of an exact dilation, should be 0): ≤ 1e−8 for M ≤ 1836, 9e−11 at 4000, 8e−7 at 10⁴. Dense windows with off-optimal fixed α have larger N−4 changes at their edges (up to 1.2e−5 for χ_F at M = 0.17 in w1) but agree with the optimized grid to ≤ 5e−6 relative.

**Rigorous vs fitted.** Rigorous: analyticity (hence no sharp features); the 1/32 limit (leading order); the small-M limits χ_F → ‖Dψ_∞‖²/(M+1)², χ_inv, χ_sep → 0. Computed to the stated basis error: all χ values and feature locations. Fitted: c₁, c₂ of the large-M expansion (window-dependent at the 10% level), and the small-M exponents.

**Interpretation.** The fidelity susceptibilities, including the two that remove trivial dilations, are smooth crossover functions. The removal of the global dilation and of the separate R/ρ dilations moves the main change of χ to M ≈ 1.4 (steepest rise of χ_inv and χ_sep) and leaves a broad, weak χ_inv maximum near M ≈ 6.9, between the L = 1 (4.21) and v = 1 (9.53) binding masses. No susceptibility singles out any of the ground-state crossover markers. The only near coincidence is the χ_sep shoulder at 0.244 with the Lindemann maximum (0.2–0.3), and both are broad. The fidelity therefore adds another descriptor-dependent crossover location (≈ 1.4) to the list rather than a common emergence point.


## Conventions

* Variable **ln M**. F(M, M e^h) = |⟨ψ(M)|ψ(M e^h)⟩| = 1 − ½ χ h² + O(h³), and
  χ = ⟨ψ̇|ψ̇⟩ − |⟨ψ|ψ̇⟩|², ψ̇ = dψ/d ln M (equivalently χ = −∂²_h ln F at h = 0 for the symmetrized F). Same as `solvers/harmonic_model` README §4, so the numbers compare directly with 3/64, 1/32, 3/16 there. To convert to the variable M divide by M²; to κ = M^(−1/4) multiply by 16/κ².
* **Fixed coordinates**: overlaps in the Jacobi coordinates R = x₁ − x₀ (X–X) and ρ = x₂ − (x₀+x₁)/2 (electron from the X–X midpoint). Both are mass independent, so this is the same as fixed interparticle distances (the variables of the solver). Coordinates measured from the centre of mass would be M dependent and are not used.
* **χ_F** (plain), **χ_inv** (global dilation removed), **χ_sep** (separate dilations of R and of ρ removed), defined by maximizing the overlap over the unitary dilations, to second order in h:
  * U_s ψ(R, ρ) = s³ ψ(sR, sρ), generator D = x·∇ + 3 (6 internal dimensions): χ_inv = χ_F − ⟨ψ̇|Dψ⟩²/‖Dψ‖² (real states, ⟨ψ|Dψ⟩ = 0).
  * U_{s,t} ψ = s^{3/2} t^{3/2} ψ(sR, tρ), generators D_R = R·∇_R + 3/2 and D_ρ = ρ·∇_ρ + 3/2 (D = D_R + D_ρ): χ_sep = χ_F − vᵀG⁻¹v with v = (⟨ψ̇|Dψ⟩, ⟨ψ̇|D_Rψ⟩) and G the Gram matrix of (Dψ, D_Rψ) (same span as {D_Rψ, D_ρψ}). This is exactly the harmonic_model definition ("separate-dilation", maximized over two independent scale factors of R and ρ).

## Method

**Wave functions and ψ̇.** For every mass the basis of the mass_scans run is rebuilt from `paper_data/mass_scans/data/L0_M*.npz` (perimetric Laguerre, exchange-symmetric 'S', the same scales α and largest N: N = 28 for M < 150, 40 for 207–1836, 44 at 4000, 48 at 10⁴ with weights (1,1,2)). The ground-state energies reproduce the stored ones to the last digit. In a fixed basis H(M) = (T₀+T₁)/M + T₂ + V, so ∂H/∂ln M = −(T₀+T₁)/M exactly, and ψ̇ (coefficients x) is the solution of the bordered system

    [[H − E S, S c], [cᵀS, 0]] [x; λ] = [−(H′ − E′S) c; 0]

(first-order perturbation theory summed exactly in the basis; xᵀSc = 0). χ_F = xᵀSx. This is the derivative of the variational state at fixed basis, which converges to the exact ψ̇ as N grows; the basis is nested (N − 4 is a sub-block with the same α), which gives the convergence check at no cost.

**Dilation generators in perimetric coordinates.** With s_k = r_kj + r_kl − r_jl, r_ij = (s_i + s_j)/2 and the volume element w = r₀₁r₀₂r₁₂/4 ds:
* D = s₀∂₀ + s₁∂₁ + s₂∂₂ + 3. A common rescaling of R and ρ rescales all three distances, hence all s_k (Euler operator). On a basis function it is exact: U_s f_n(s; α) = f_n(s; sα), so overlaps with dilated states are available exactly (`fidcore.dilated_overlap`).
* D_R = R·∇_R|_ρ + 3/2. Since r₀₂² − r₁₂² = 2ρ·R, one finds R·∇_R r₀₁ = r₀₁, R·∇_R r₀₂ = (r₀₁² + r₀₂² − r₁₂²)/(4r₀₂), R·∇_R r₁₂ = (r₀₁² − r₀₂² + r₁₂²)/(4r₁₂), so

      D_R = r₀₁∂_{r01} + (r₀₁²+r₀₂²−r₁₂²)/(4r₀₂) ∂_{r02} + (r₀₁²−r₀₂²+r₁₂²)/(4r₁₂) ∂_{r12} + 3/2,
      ∂_{r_ij} = ∂_{s_i} + ∂_{s_j} − ∂_{s_k}.

  Separate dilation of R and ρ is **not** a rescaling of the perimetric coordinates (it changes the triangle's shape), so it has no closed form in the basis; the generator approach is exact anyway because only the scalars ⟨ψ̇|D_Rψ⟩, ⟨Dψ|D_Rψ⟩, ‖D_Rψ‖² are needed. The bounded coefficients (e.g. (r₀₁²−r₁₂²)/r₀₂, |r₀₁−r₁₂| ≤ r₀₂) make D_Rψ smooth except at the coalescence corners.
* **Integration.** All integrands except two are polynomial × exp(−Σα_k s_k): product Gauss–Laguerre with n_max+8 nodes per coordinate is exact (checked: n_max+16 changes nothing at 1e−15). The two singular pieces of ‖D_Rψ‖², ∫(r₀₁r₁₂/4r₀₂)B̃² and ∫(r₀₁r₀₂/4r₁₂)C̃², are done in polar coordinates (s₀, s₂) = (tu, t(1−u)) [resp. (s₁, s₂)], where 1/r₀₂ = 2/t cancels the Jacobian t: exact Gauss–Laguerre in t for each u, composite Gauss–Legendre in u graded to both ends (14 panels per half × 16 nodes; 20 × 24 or more t nodes change the result by < 1e−14).

## Validation (`validation.txt`; `validate.py M`, `validate_DR.py M`, `limit_small_M.py`)

At M = 1 (N = 28):
* χ_F, perturbative 0.469825574773; finite differences of the fixed-basis overlaps F(M, Me^{±h}) (averaged ±h), h = 0.02: 0.469806434, h = 0.01: 0.469820790, Richardson 0.469825574844 (agreement 1.5e−10 relative). The FD check is stored for every grid mass (`chiF_fd_rich`).
* Quadrature norm 1 − 1e−14; ‖ψ̇‖² by quadrature equals xᵀSx to 12 digits; anti-Hermiticity ⟨ψ|Dψ⟩, ⟨ψ|D_Rψ⟩ ≈ 1e−15 (a non-trivial check of the D_R coefficients and of the 3/2); the two singular pieces are equal to 12 digits (exchange symmetry).
* ‖Dψ‖² = 1.4670745112 against 1.4670745108 from exact dilated overlaps (σ = 0.01, 0.005, Richardson).
* **Coulomb mass-scaling identity**: multiplying all three masses by e^h is an exact dilation, ψ(h) = U_{e^h}ψ. The fixed-basis derivative for that perturbation gives ‖ψ̇‖² = 1.46707449 vs ‖Dψ‖² = 1.46707451, ⟨ψ̇|Dψ⟩/‖Dψ‖² = 1 − 9e−9, and χ_inv = χ_sep = 5.5e−9 (exact 0). Stored for every mass as `scal_chi_inv` (a direct measure of how well the basis represents a dilation).
* **Finite-separation dilation-invariant fidelity**: max over σ of |⟨ψ(M)|U_σψ(Me^{±h})⟩| (exact mixed-scale overlaps, 1D maximization): h = 0.02 → 0.0099818642, h = 0.01 → 0.0099815266, h² extrapolation 0.0099814141, equal to the projected χ_inv = 0.0099814141.
* **Independent brute-force check of D_R**: direct quadrature in Jacobi coordinates (R, ρ, θ) with D_Rψ from a 4th-order finite difference in R at fixed ρ, θ: ‖D_Rψ‖², ⟨ψ̇|D_Rψ⟩, ⟨Dψ|D_Rψ⟩ agree with the perimetric result to 1.3e−6, 1.2e−6, 7.5e−7 (the accuracy of the brute-force quadrature).
* **Small-M limit**: as M → 0 the state tends to a dilated infinite-mass H⁻-like state (two light X around a static e), ψ(M) ≈ U_μ ψ_∞, μ = M/(M+1), so χ_F → ‖Dψ_∞‖²/(M+1)² and χ_inv, χ_sep → 0. ‖Dψ_∞‖² = 1.3675464 (masses (1,1,∞), E = −0.5277510165, N = 24 and 28 agree to 4e−9).

## Files

| file | content |
|---|---|
| `fidcore.py` | engine: basis rebuild, ground state, bordered ψ̇ solve, FD check, perimetric quadrature of D, D_R scalars, exact dilated overlaps |
| `run.py` | resumable driver: `grid` (mass_scans masses → `results.csv`), `dense w1..w4` (fixed-α windows w1–w7 → `dense.csv`; w1 0.17–0.40 (α from M = 0.3), w2 0.33–0.66 (0.5), w3 0.5–0.95 (0.7), w4 0.8–1.35 (1.0), w7 1.2–2.6 (1.5), w5 3–12.5 (6), w6 11.5–50 (20)) |
| `validate.py`, `validate_DR.py`, `limit_small_M.py` | checks above, appended to `validation.txt` |
| `features.py` | locations of extrema and inflections with N/N−4 and interpolation spread → `features.txt` |
| `analyze.py` | `fits.txt`, `fig1_chi_vs_M.png`, `fig2_dense_windows.png`, `fig3_largeM_convergence.png` |
| `results.csv`, `dense.csv` | one row per mass: χ_F (perturbative and FD), χ_inv, χ_sep, the scalars ‖Dψ‖², ‖D_Rψ‖², ⟨Dψ|D_Rψ⟩, ⟨ψ̇|Dψ⟩, ⟨ψ̇|D_Rψ⟩ (all per ‖ψ‖²), checks, the same three χ at N − 4, the uniform-scaling test |
| `data/dpsi_<key>.npz` | c and x (ψ and ψ̇ coefficients), α, N, weights: reuse without recomputation |
| `log_*.txt` | run logs |
