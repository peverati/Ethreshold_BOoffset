# Why BO critical masses are offset by an almost constant amount

*Last updated: 2026-11-14 17:08*

**Update 2026-11-14:** the mass functions m_v(R), m_r(R) are now computed from first principles in `massfunc/` and replace the fitted effective masses: see `massfunc/README.md` and `massfunc/analysis.txt`.

Written 2026-11-14. System: {X+, X+, e-}, nuclear mass M in electron masses. Question: why the BO (and BO plus adiabatic) critical masses exceed the exact ones by almost the same amount for every vibrational level (1.53 to 1.56 for BO, 0.50 to 0.54 for BO+ad), and why the offset is smaller for L = 1 and L = 2. Exact M_c come from `paper1/critical_masses/` and `paper1/rotational/`. Everything else was computed here with `solvers/bo_h2plus` tables (U(R), D(R)). No new three-body runs were made.

## 1. Main result (derived)

**The critical mass is an eigenvalue of a linear problem.** At threshold the BO radial equation becomes a zero-energy equation. Two choices make it simple:

* scale R by the electron reduced mass, r = mu_e R with mu_e = 2M/(2M+1);
* multiply the equation by M/mu_e².

The result is, exactly,

    [-d²/dr² + L(L+1)/r² + c ΔD(r)] u = λ (−ΔU(r)) u,     λ = M/mu_e = M + 1/2,

where:

* ΔU = U(r) + 1/2 is the unit-mass BO potential counted from its asymptote;
* ΔD = D(r) − 1/4 is the adiabatic correction counted from its asymptote;
* c = 0 gives BO and c = 1 gives BO+ad.

A critical mass is an eigenvalue: **M_c = λ_c − 1/2.** The coupling M enters linearly, so a correction that scales like the nuclear kinetic energy (order 1/M in the energy) becomes an **M-independent perturbation** of this eigenproblem. Such a correction therefore shifts M_c by an amount that does not depend on M.

* **Check.** `coupling.py` (shooting, with the exact −C4/r⁴ zero-energy tail used for matching) reproduces every critical mass of `bo_h2plus` to ≤ 1e−9: 0.9649154241, 11.0814565935, 31.1677742733, 60.8466743844 (BO); 10.0631790317, 30.1455551746, 59.8156816748 (BO+ad); L = 1: 5.5374983483 and 4.5705652790; L = 2: 12.7397228367 and 11.7786163976.
* **Why the naive route fails.** Linearizing the binding energy, ΔM ≈ ΔE/(dE_b/dM), does not work at an L = 0 threshold. There E_b ∝ (M − M_c)², so dE_b/dM = 0 at M_c. The coupling-constant formulation is the correct linear version of the same idea.
* **First-order perturbation theory for λ.** Write A = −d² + L(L+1)/r² and B = −ΔU. A perturbation δA shifts the eigenvalue by δλ = ⟨u|δA|u⟩/⟨u|B|u⟩ (u is the zero-energy solution). This is ordinary Rayleigh–Schrödinger theory for the generalized eigenproblem; the integrals converge because B ∝ r⁻⁴.

**Why the shift is almost independent of v.** In the semiclassical (WKB) limit, u² ≈ 1/p with p = sqrt(λ(−ΔU)). Every zero-energy level then samples the same weight, so δλ tends to a v-independent limit. The remaining v-dependence is a quantum correction that decreases as v grows. Example (adiabatic shift): 1.016, 1.021, 1.030 for v = 1, 2, 3, with WKB limit 1.058 (`wkb_limits.txt`). The weight does **not** sit in the halo tail: more than 99% of it lies at r < 8 bohr (`weights_cumulative.csv`). So the offset measures the corrections inside the molecular region, not at dissociation.

## 2. Decomposition of the offset

**(a) BO → BO+ad (derived, no fit).** In first order,

    M_c^BO − M_c^BOad ≈ ⟨1/4 − D⟩_u / ⟨−ΔU⟩_u,

the adiabatic correction relative to its asymptote, divided by the well depth, both averaged over the zero-energy solution. It gives 1.0163, 1.0212, 1.0304 (L = 0, v = 1, 2, 3), 0.9655 (L = 1) and 0.9601 (L = 2). The exact differences (both models solved) are 1.0183, 1.0222, 1.0310, 0.9669 and 0.9611, so first order is good to 0.2%. The value is about 1 because D(∞) − D(R) (about 0.075 at R_e) and the well depth (about 0.10) are comparable in the region the state samples.

**(b) BO+ad → exact: the effective vibrational mass (derived structure, one input fitted).** The leading nonadiabatic effect is an R-dependent kinetic mass: M → M + m_v(R) for radial motion and M → M + m_r(R) for rotation. A constant m shifts λ by

    δλ = −[m_v ∫u'² + m_r ∫L(L+1)u²/r²] / (λ ∫(−ΔU)u²) = −m_v f_rad − m_r f_cen.

For L = 0, f_rad = 1 + |δλ_ad|/λ = 1.097, 1.033, 1.017 (the 1/λ part arises because the adiabatic term is not rescaled by the mass).

* **Asymptotic value (derived).** At R → ∞ the relative motion is that of X+ against the atom (Xe), with reduced mass μ = M(M+1)/(2M+1). Then 2μ = M + 1/2 − O(1/M), so m_v(∞) = m_r(∞) = 1/2.
* **Values near R_e (from data, this folder).**
  * m_v(R_e) = 0.3465 ± 0.0001. It is the κ⁶ coefficient of the ground-state energy error (see `kappa6.md`).
  * m_r(R_e) = 0.0754 ± 0.001. It comes from the exact and BO+ad J = 1, 2 spacings, −M(d − d_BOad)/d_BOad → 0.0754 as M → ∞ (`rot_mass_data.csv`).

## 3. Test table

(`offset_table.csv`.) "pred nonad" uses m_v = 1/2 (no fit) and, for L > 0, m_r = m_r(R_e) = 0.0754. "pred BO" = exact (BO − BO+ad) + pred nonad.

| state | M_c exact | offset BO | offset BO+ad | pred BO | pred nonad (BO+ad) | residual (BO+ad) | m_eff from data |
|---|---|---|---|---|---|---|---|
| L=0, v=1 | 9.525634 | 1.5558 | 0.5375 | 1.5666 | 0.5483 | −0.0108 | m_v = 0.4902 |
| L=0, v=2 | 29.63939 | 1.5284 | 0.5062 | 1.5389 | 0.5167 | −0.0105 | m_v = 0.4898 |
| L=0, v=3 | 59.31866 | 1.5280 | 0.4970 | 1.5395 | 0.5086 | −0.0115 | m_v = 0.4887 |
| L=1 | 4.2082827 | 1.3292 | 0.3623 | 1.2334 | 0.2665 | +0.0958 | ⟨m_r⟩ = 0.205 (with m_v = 0.4895) |
| L=2 | 11.4894217 | 1.2503 | 0.2892 | 1.1463 | 0.1852 | +0.1040 | ⟨m_r⟩ = 0.203 (with m_v = 0.4895) |

**Reading of the table.**

* **L = 0.** Using m_v = 1/2 with no fit gives the BO+ad offset to 0.011, and the BO offset to 0.011 (0.7%). The residual is the same for all three levels. Expressed as an effective vibrational mass, the data give m_v = 0.4895 ± 0.0008, constant over M_c = 10 to 60.
  * Constancy implies that any O(1/M²) nonadiabatic potential term (which would enter as b/M_c) is small: |b| ≲ 0.02.
  * The total BO offset is roughly constant because two near-constant pieces add: 1.02 (adiabatic) + 0.49 × (1 + 1.02/λ) (mass). Their small v-trends have opposite signs (1.018 → 1.031 and 0.538 → 0.497) and largely cancel.
* **L > 0.** The offsets are smaller because a large fraction of the kinetic weight is centrifugal (f_cen = 0.77 for L = 1, 0.83 for L = 2), and the rotational mass correction is much smaller than the vibrational one near R_e.
  * m_v = m_r = 1/2 overestimates the offsets (0.595, 0.539).
  * m_r = m_r(R_e) underestimates them by 0.10 for both.
  * Fitting a single effective ⟨m_r⟩ to L = 1 (0.205) gives 0.203 for L = 2 (difference 0.0015 in m_r, i.e. 0.0012 in M_c). That consistency is a test of the structure, but **⟨m_r⟩ itself is fitted**.

## 4. What is derived, fitted, or unexplained

* **Derived (exact or first-order theory, checked numerically):**
  * the linear eigenvalue form with λ = M + 1/2;
  * why any O(1/M) correction gives an M-independent shift;
  * why it is nearly v-independent (WKB);
  * the adiabatic part of the offset (1.02 and 0.96–0.97), to 0.2%;
  * the asymptotic masses m_v(∞) = m_r(∞) = 1/2;
  * the 1/λ factor f_rad.
* **Taken from data (fits):**
  * m_v(R_e) = 0.3465 and m_r(R_e) = 0.0754 (large-M fits);
  * the effective values ⟨m_v⟩ = 0.4895 and ⟨m_r⟩ = 0.204 (from the offsets themselves).
* **Unexplained / open:**
  1. ⟨m_v⟩ = 0.49 is a weighted mean of the unknown function m_v(R) over r ≈ 1.5–8 bohr. Yet m_v(R_e) = 0.35, and 13–36% of the weight lies at r < 3. So m_v(R) must rise quickly beyond R_e and probably exceed 1/2 at intermediate R (a hump), or other nonadiabatic terms contribute a near-constant shift. This is our inference and untested. The decisive test is to compute the H2+ nonadiabatic mass functions m_v(R) and m_r(R) (second-order perturbation theory; Dalgarno–Lewis response of ∂_Rψ, gerade only) and fold them with `weights_cumulative.csv`. This was not done, for time reasons.
  2. The L > 0 residuals (+0.10) are reproduced only with a fitted ⟨m_r⟩ ≈ 0.20, which is not computed from m_r(R). For L = 1 (M_c = 4.2), 1/M_c corrections may be sizeable and are absorbed into that fit.
  3. Literature values of m_v(R), m_r(R) for H2+ (Bunker–Moss type vibrational and rotational masses; Kutzelnigg, "which masses are vibrating or rotating", Mol. Phys. 2007; Jaquet and Kutzelnigg 2008) were **not checked**. Our recollection that the "vibrational mass is closer to the atomic mass and the rotational mass closer to the nuclear mass" is **from memory, unverified**. It agrees qualitatively with m_v ≫ m_r found here.
* **Convention note.** The BO offset depends on the BO convention. bo_h2plus BO already uses mu_e scaling. A clamped-nuclei BO with electron mass 1 and nuclear kinetic −(1/M)∇² has λ = M, so its critical masses are 1/2 larger (offset ≈ 2.03 for L = 0). The BO+ad offset does not depend on this choice to the order considered.

## Files

* `coupling.py`: zero-energy coupling eigenproblem (BO, BO+ad) and first-order integrals. Output: `coupling_results.csv`.
* `weights.py`: radial distribution of the weights and the WKB limits. Output: `weights_cumulative.csv`, `wkb_limits.txt`.
* `offset_table.py`: the test table. Output: `offset_table.csv`.
* `rot_mass_data.csv`: m_r(M) from exact and BO+ad J-spacings (exact from `mass_scans/exact_rot.csv`, BO+ad from `mass_scans/bo.csv`).
* `kappa6.md`, `kappa6_data.csv`, `fit_kappa6.py`, `kappa6_fits.csv`: Task B (the κ⁶ coefficient).
* `PROGRESS.md`: step log.

Run from the repo root with `python3 paper1/bo_offset/<script>.py`. Total runtime is under 2 minutes.
