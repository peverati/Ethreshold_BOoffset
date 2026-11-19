# The κ⁶ coefficient of E_exact − E_BO+ad

*Last updated: 2026-11-14 15:42*

Written 2026-11-14. Ground state (L = 0, v = 0) of {X+, X+, e-}, κ = M^(−1/4). E_BO+ad is the solvers/bo_h2plus convention:

* nuclear kinetic energy −(1/M)d²/dR², i.e. the **nuclear** reduced mass M/2;
* potential mu_e U(mu_e R) + mu_e² D(mu_e R)/M, with mu_e = 2M/(2M+1).

## Data

* Exact energies: `paper1/mass_scans/energies.csv` (E0_extrap).
* E_BO+ad: same file (E_BOad), checked against `mass_scans/bo.csv`.
* M = 10⁴: E = −0.60032852464759 (±8e−13) from `paper1/convergence_large_M/results.csv`. It coincides with the mass_scans value (−0.600328524647589) to 1e−15, so the fits do not change.
* Fit table: `kappa6_data.csv` (M ≥ 20). The mass_scans CSVs were not modified.

## Coefficient (numerical, solid)

Fixed-exponent fits y = M^(3/2)(E − E_BO+ad) = c6 + c8 κ² (+ c10 κ⁴), weighted by the exact-energy error, over several windows (`fit_kappa6.py` → `kappa6_fits.csv`):

| fit | windows (M_min) | c6 | c8 |
|---|---|---|---|
| 3 terms | 100, 206.8, 500, 1000 | −0.039278, −0.039303, −0.039311, −0.039315 | −0.0760 … −0.0740 |
| 2 terms | 500, 1000, 1836 | −0.039335, −0.039327, −0.039320 | −0.0728 … −0.0734 |

**c6 = −0.03931 ± 0.00003**. The error is the spread over the 3-term windows with M_min ≥ 206.8 and the 2-term windows with M_min ≥ 1000. **c8 ≈ −0.074 ± 0.002** (less certain: it depends on the fit order). The earlier two-term value in mass_scans/fits.txt (−0.039339, M ≥ 500) is consistent.

## Physical origin (derived structure, number from the fit)

* **Order counting.** The order-κ⁶ = M^(−3/2) term of the nonadiabatic correction is the **R-dependent vibrational mass**. Second-order perturbation in the nonadiabatic coupling gives an operator (1/M²) ∂_R a(R) ∂_R; acting on a vibrational function with ⟨p²⟩ ∝ M^(1/2), it gives order M^(−3/2). Nonadiabatic *potential* terms are O(M^(−2)) = κ⁸.
* **Mass form.** Writing the correction as M → M + m_v(R) in the radial kinetic energy, first order gives E − E_BO+ad ≈ −(m_v/M)⟨T_rad⟩. In the harmonic limit ⟨T_rad⟩ = ω/4 and √M ω/4 = sqrt(2k)/4 = 0.113452, with k = U''(R_e) = 0.1029705 (unit mass). So

      c6 = −m_v(R_e) · sqrt(2k)/4   ⇒   m_v(R_e) = 0.3465 ± 0.0003 electron masses.

* **Independent check.** m_v_eff(M) = −M(E − E_BO+ad)/⟨T_rad⟩_BO+ad, with the BO+ad ⟨T_rad⟩ per mass, extrapolated in κ² (3-term fits) gives **0.34654 ± 0.0001**, the same value.
* **Meaning.** The nuclei vibrate as if each carried about 0.35/2 ≈ 0.17 of an electron mass extra (in the convention of the nuclear reduced mass M/2). The separated-atom limit is m_v(∞) = 1/2: the reduced mass of X+ and (Xe) is M(M+1)/(2M+1) ≈ (M + 1/2)/2.
* **Rotational analogue.** m_r(R_e) = 0.0754 (from J = 1, 2 spacings; README §2), much smaller than m_v.
* **Link to Task A.** The same m_v, averaged over the zero-energy critical states, gives ⟨m_v⟩ = 0.4895, which is the BO+ad critical-mass offset (README).

## Convention dependence (derived)

* **Potential convention: no effect at this order.** The two are:
  * mu_e-scaled BO + adiabatic (used here);
  * unscaled U(R) plus an adiabatic correction that includes the electron's mass-polarization/reduced-mass term.

  They agree to O(1/M) in the potential, because mu_e U(mu_e R) = U + T_el(R)/(2M) + O(M⁻²) by the virial relation. They differ at O(M⁻²), which changes the ground-state energy only at O(M^(−5/2)). **c6 does not depend on this choice.**
* **Kinetic-mass convention: shifts c6 directly.** Using reduced mass m_N/2 with m_N = M + δ shifts m_v by −δ, so c6 → c6 + δ · 0.113452:
  * nuclear mass (δ = 0, used here): c6 = −0.03931;
  * "nucleus plus half an electron" (δ = 1/2): c6 = +0.01742;
  * atomic mass (δ = 1): c6 = +0.07414.

  Any quoted κ⁶ coefficient must therefore state the nuclear-mass convention.

## Literature comparison

Not done. We have not checked any published H2+ nonadiabatic vibrational mass function (Bunker–Moss type, or Kutzelnigg 2007 / Jaquet and Kutzelnigg 2008). Any number I could quote for them would be **from memory and unverified**, so none is quoted. A consistency check that can be made:

* At H2+ (M = 1836.15) the full E − E_BO+ad is −5.215e−7 Eh = −0.1145 cm⁻¹.
* The κ⁶ term alone gives c6 M^(−3/2) = −5.00e−7; c8 κ⁸ adds −2.2e−8.

This is a decomposition of our own data, not an external check.

## What stays open

* The functions m_v(R) and m_r(R) themselves (needed for Task A, item 1 of README §4) are not computed.
* c8 mixes the anharmonic and vibrational-averaging part of m_v(R) with the nonadiabatic potential term. It is not interpreted.
