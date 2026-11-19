# Comparison with King, Longford & Cox, JCP 139, 224306 (2013)

*Last updated: 2026-11-19 13:35*

Written 2026-11-13. This folder checks our exact three-body solver (`solvers/threebody`, perimetric-Laguerre basis, used unchanged) against all ground-state energies of the symmetric family {m1, m1, m3} in King et al. (Table III: 42 finite-mass rows plus infinite-mass H⁻, 43 in all; Table II: 2856-term benchmarks and their "literature" column). It also locates the minimum of their marker g = (E0 − E_th)/E_th without fitting.

Source: the paper itself (CC BY, Sussex figshare) and its Publisher's Note. Table III has 43 rows (42 finite masses plus ∞H⁻).

## Files

| file | content |
|---|---|
| `king2013_table3.csv` | Table III transcribed as printed (energy with its uncertain digits in parentheses, m1/m3, a3, g). The `note` column records the ambiguities. |
| `king2013_table2.csv` | Table II: all 1078-term variants, the 2856-term values and the literature column, with notes |
| `common.py` | CODATA-2012 masses (their Table I), parser for "−0.599506(9)", g, a3 |
| `compute.py` | stages `table`, `fine`, `minimum`, `extra`, `lit` (restartable; each appends rows) |
| `analyze.py` | builds the comparison tables, fits, minimum and figure → `results.txt` |
| `ours_table.csv` | our E0 for the 38 distinct mass ratios: alphas, nested N sequence `E_seq`, extrapolation, error, flag, virial |
| `ours_extra.csv` | 7 extra points a3 = 0.17 … 0.23 (the extra points King added for their fit; their exact a3 values are not printed) |
| `ours_lit.csv` | Table II systems at the masses behind their "literature" column |
| `fine_grid.csv` | 201 points, M = 0.300 … 0.800, step 0.0025: E, the exact dE/dM, g and dg/dM at N = 28 and 32 |
| `minimum.csv` | root of dg/dM = 0 for four scale choices and N = 20 … 32 |
| `comparison_table3.csv`, `comparison_table2.csv` | ours − theirs, row by row |
| `results.txt` | all printed results of `analyze.py` |
| `fig_g_comparison.png` | (a) g(M) ours and theirs; (b) the minimum region, with the dense grid, their points, our minimum, their 0.491 ± 0.003 band and the Eq. 9 fit minimum; (c) relative difference \|ours − theirs\|, their last certain digit, and our error estimate vs M |
| `run_log.txt`, `log_*.txt` | logs |

Regenerate with `export THREEBODY_BLAS_THREADS=3; python compute.py table; python compute.py fine; python compute.py minimum; python compute.py extra; python compute.py lit; python analyze.py`, using `python3`. Wall time was about 1.5 h with two processes on a loaded machine.

## Method

**Systems and units.** Every system in both tables has two identical particles of mass m1 and one particle of opposite charge with mass m3. Charge inversion does not change the energy. That includes H⁻ = (e,e,p) and Mu⁻ = (e,e,μ), which are the {M, M, 1} family with M = m1/m3 < 1, so no separate "inverse arrangement" code was needed. Each system is computed as one dimensionless problem with the smaller of m1 and m3 set to 1: masses (M, M, 1) for M ≥ 1, (1, 1, 1/M) for M < 1, and (1, 1, ∞) for ∞H⁻. The result is then scaled by that mass in m_e to give King's electron atomic units. Their muonic-unit entries in Table II (μH⁻, μD₂⁺) are exactly our dimensionless energies, because the muon is the unit particle there. Masses are CODATA 2010 (Mohr, Taylor and Newell, published 2012; the folder names say "2012"), as in their Table I (p 1836.15267245, d 3670.4829652, t 5496.9215267, μ 206.7682843, τ 3477.15). The hypothetical rows are exactly M = 11.3, 6, 4.7, 3.7, 2.4 and 0.2. I checked this: g recomputed from their printed E with these M matches their g column to ≤ 2e−8.

**Basis and convergence.** The settings are those of `paper_data/mass_scans`, with one extra N for light systems. Scales (α_X, α_e) are optimized by Nelder–Mead at N_opt and then **fixed**. The sequences are nested, so the energies decrease monotonically:
- M < 150: total degree, N_opt 12, N = 16, 20, 24, 28, 32 (size 3417)
- 150 ≤ M < 3000: weights (1,1,2), N_opt 20, N = 28 … 40
- 3000 ≤ M < 5000: N = 32 … 44
- M ≥ 5000: N = 36 … 48 (size 5525)

Extrapolation uses the package's geometric-tail rule. Its error is max(tail, change of the extrapolated value), which is a heuristic and not a bound. Errors in the dimensionless energy are ≤ 2.7e−11 for every mass, ≤ 1e−12 for M ≥ 3, and 4e−12 at the heaviest masses. The light (e,e,X) systems converge most slowly, with a last increment of about 4e−11 at N = 32. In King's units the error is multiplied by the unit mass, for example 5e−8 hartree for (t,t,t) ≈ −1440 hartree. That is still 100 times below their last certain digit.

**Reproducibility check.** Where `paper_data/mass_scans` uses the same M, the two runs agree within the error bars: M = 1 differs by 6e−12, M = 6 by 1.7e−13, and M = 0.2 by 1.8e−11. The muon and proton masses differ between the two folders (CODATA 2012 here, CODATA 2018 there), so no other values were reused.

## Results

### Table III (42 finite-mass rows plus ∞H⁻)

**All 43 rows agree with King et al. within one unit of their last certain digit.** In units of that digit, |ours − theirs| ≤ 0.16 for 42 rows and 0.64 for (6,6,1). In 42 of 43 rows, truncating our value reproduces their certain digits exactly. The exception is (6,6,1): their −0.466087497252(4) against our −0.466087497253040 ± 2e−13, which is consistent with rounding at their uncertain digit. Their stated precision is conservative everywhere. Our g agrees with their g column to ≤ 1.4e−8, and to about 1e−10 for most rows.

Some sample rows (their units):

| system | M | King | ours ± err |
|---|---|---|---|
| (e,e,p) H⁻ | 5.4462e−4 | −0.527445881(0) | −0.527445881101 ± 2e−11 |
| ∞H⁻ | 0 | −0.52775101(6) | −0.527751016536 ± 2e−11 |
| (e,e,e) Ps⁻ | 1 | −0.262005070(2) | −0.262005070230 ± 9e−12 |
| (p,p,d) | 0.50025 | −639.557763(2) | −639.557763236 ± 4e−8 |
| (p,p,μ) | 8.8802 | −102.22351295(5) | −102.223512955560 ± 2e−11 |
| (μ,μ,e) Mu₂⁺ | 206.77 | −0.585126098252(3) | −0.585126098252403 ± 4e−15 |
| (p,p,e) H₂⁺ | 1836.15 | −0.59713906(3) | −0.597139063079291 ± 2e−13 |
| (t,t,e) | 5496.9 | −0.599506(9) | −0.599506910096 ± 1e−12 |

### Table II (their 2856-term benchmarks)

| system | ours (CODATA 2012) | King 2856 | ours − King |
|---|---|---|---|
| H⁻ | −0.5274458811007 ± 2e−11 | −0.527445881063 | −3.8e−11 |
| Mu⁻ | −0.5250548065238 ± 2e−11 | −0.525054806510 | −1.4e−11 |
| μH⁻ (μ.a.u.) | −0.4718663360028 ± 3e−11 | −0.471866336008 | +5e−12 |
| Ps⁻ | −0.2620050702296 ± 9e−12 | −0.262005070191 | −3.9e−11 |
| μD₂⁺ (μ.a.u.) | −0.5311111302292 ± 1e−14 | −0.531111130229 | −2e−13 |
| Mu₂⁺ | −0.5851260982524 ± 4e−15 | −0.585126098252 | −4e−13 |
| H₂⁺ | −0.5971390630793 ± 2e−13 | −0.597139063079 | −3e−13 |
| D₂⁺ | −0.5987887843048 ± 4e−12 | −0.598788784032 | **−2.7e−10** |

The molecular-type systems Mu₂⁺, H₂⁺ and μD₂⁺ agree to the 12 digits they print. For the atomic-type H⁻, Ps⁻ and Mu⁻ and for D₂⁺, their 2856-term value lies above ours, by 1–4e−11 for the first three and 2.7e−10 for D₂⁺. Their method is variational and so is ours, so the lower value is the better upper bound. These gaps are therefore incomplete convergence on their side, which the paper itself says is slower for molecular systems. It does not affect any Table III entry, because Table III quotes fewer digits.

**Their "literature" column is explained by the older mass constants.** I recomputed each system with the masses listed in the caption of their Table I for the cited works (`ours_lit.csv`). Ours − literature: H⁻ +3.5e−12, Mu⁻ +8e−12, Ps⁻ +2.4e−12 (mass-independent), Mu₂⁺ −3e−13, H₂⁺ −1.7e−13 (m_p = 1836.152701), D₂⁺ −1.3e−13 (m_d = 3670.483014), μD₂⁺ +4.3e−11. Two loose ends remain:
- **μH⁻ is off by 6.4e−9.** The literature value is printed with only 10 digits and without its minus sign, and I may not have the mass set of their Ref. 39 right. This one is unresolved.
- **μD₂⁺ is off by 4e−11.** Frolov's deuteron mass is taken here to equal theirs, which may not be exact.

**A side result for the solver.** `solvers/threebody/README.md` lists an "open discrepancy at the 1e−10 level" for H₂⁺ against −0.597139063123 (a value quoted from memory). That number is the H₂⁺ energy for m_p = 1836.152701. At that mass our solver gives −0.5971390631232, so **the discrepancy is just the proton mass and is closed.** This is the earlier README's mistaken estimate of the mass sensitivity. I did not edit the solver README; that is for you to decide.

### The minimum of g(M)

g = (E0 − E_th)/E_th with E_th = −M/(2(M+1)), in their units (odd particle mass 1). In King's coordinate, a3 = M/(M+2).

**Fit-free (rigorous apart from the basis):** **M_min = 0.4924222 ± 0.0000001**, which is a3 = 0.1975678, with **g_min = 0.0451139164** and curvature d²g/dM² = 0.0380. Two independent routes give this:
1. **Brent root of dg/dM = 0.** dg/dM uses dE/dM = −(⟨T₀⟩ + ⟨T₁⟩)/M², the exact Hellmann–Feynman derivative of the variational eigenvalue at fixed basis. Results (`minimum.csv`):
   - N = 20, 24, 28, 32 with scales optimized at M = 0.49: 0.49242282, 0.49242242, 0.49242231, 0.49242227. The sequence converges geometrically.
   - Four different scale sets at N = 32 spread over 0.49242226 … 0.49242232 (6e−8).
2. **Dense grid** (201 points, step 0.0025, scales re-optimized every 0.05):
   - The grid minimum is at 0.4925 ± 0.00125 (half step).
   - A 3-point parabola gives 0.492426.
   - The sign change of the exact dg/dM falls at 0.492423.
   - N = 28 and N = 32 give identical results, with E(32) − E(28) ≤ 2.8e−11 over the grid.

The basis error in M_min (≤ 1e−7) follows from the energy errors (~1e−11) and the curvature. The quoted ±1e−7 covers the N trend and the scale spread.

**In King's fit form.**
- **Their printed coefficients cannot be used.** Their Eq. 9 is g(a3) = Σ_{i=0}^{10} c_i (1 − a3)^{i/2}. At a3 = 1 every term except c0 vanishes, so g(1) = c0 = 0.5886 against data ≈ 0.20. With all coefficients positive as printed, g(0.2) = 261.
- **Refitting their form reproduces their fit statistics.** Least squares on the 37 distinct finite-mass points plus ∞H⁻ gives R² = 0.99999998, the same as theirs. With their printed g the minimum is M = 0.49025 (a3 = 0.19687). With our exact g it is 0.49025, and 0.49040 with our 7 extra a3 points added.
- **The printed values look like magnitudes with lost minus signs.** Our refitted coefficients alternate in sign: +0.2053, −0.326, +0.369, −1.88, +11.4, −43.3, +102.6, −152.2, +137.5, −69.2, +14.9. Their printed magnitudes c1 … c10 (0.327, 0.383, 2.05, 12.4, 46.9, 110, 163, 146, 73.3, 15.7) are within about 10% of these, but c0 does not match (0.5886 against 0.2053). This is consistent with lost minus signs plus a misprinted c0, but I cannot prove what they used.
- **The fit result depends on the order, not on the data accuracy.** With our exact data the minimum moves with the number of terms: 7 terms 0.5187, 9 terms 0.4888, 11 terms (theirs) 0.4903, 13 terms 0.49238. Their simpler Eq. 8, A − B√(1−a3) + C(1−a3)², has R² = 0.998, the same as theirs (A, B, C = 0.2010, 0.2379, 0.0915, against their printed 0.200907, 0.237411, 0.091054). Its minimum sits at M = 0.666, far off.

**Verdict.**
- King's M = 0.491 ± 0.003 (a3 = 0.197 ± 0.001) contains the exact value 0.49242, which sits about half their error bar from their centre.
- Their central value is low by 0.0014. Their own form, refitted to our exact data, is low by 0.0020 (0.4904 against 0.4924). That bias comes from the truncated √(1−a3) series, not from their energies, which are accurate enough.
- **What is rigorous:** M_min = 0.4924222(1), as the minimum of the variational g at the converged basis. **What is fitted:** every value from Eq. 8 or 9, which carries a form bias of up to ±0.03 depending on the order.

## Discrepancies and likely causes (summary)

1. **Table III typos.**
   - (τ,τ,e): m1/m3 is printed 3.477500e3. Their own Table I (3477.15) and their a3 column both imply 3477.15.
   - (4.7,4.7,1): the minus sign is missing from the energy.
   - Hypothetical rows: a3 is printed rounded (0.85, 0.70, 0.65, 0.55, 0.09). The M values are exact.
2. **Eq. 9 coefficients are misprinted.** See above.
3. **Table II convergence:** their 2856-term D₂⁺ is 2.7e−10 too high, and H⁻, Mu⁻ and Ps⁻ are 1–4e−11 too high (variational argument). Table III is unaffected.
4. **Their "literature" column** differs because of older mass constants, as they say. It is reproduced to ≤ 1e−11 except for μH⁻ (6.4e−9, unresolved) and μD₂⁺ (4e−11).
5. **"36 distinct finite-mass systems"** in their text: Table III has 37 distinct finite mass ratios, so their count presumably excludes one row. This does not matter for the fit.

Reference values not from this paper (best Ps⁻ ≈ −0.26200507023298, ∞H⁻ ≈ −0.527751016544) are from memory and only corroborate. Ours agree with them to about 1e−11, and nothing above depends on them.
