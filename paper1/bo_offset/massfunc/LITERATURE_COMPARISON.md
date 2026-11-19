# Literature comparison: nonadiabatic vibrational and rotational mass functions
*Last updated: 2026-11-15 21:25*

Written 2026-11-15. Compares our m_v(R), m_r(R) (`results_fine.csv`, `THEORY.md`) with three papers read in full:

* [BM] P. R. Bunker and R. E. Moss, Mol. Phys. **33**, 417 (1977), doi:10.1080/00268977700100351.
* [BMM] P. R. Bunker, C. J. McLarnon and R. E. Moss, Mol. Phys. **33**, 425 (1977), doi:10.1080/00268977700100361.
* [K] W. Kutzelnigg, Mol. Phys. **105**, 2627 (2007), doi:10.1080/00268970701604671.

Sections 7–9 add Moss 1996 [M96], Schwenke 2001 [S01] and Jaquet–Kutzelnigg 2008 [JK08]; Section 10 is the current novelty statement and supersedes Sections 4 and 7 on that point.

Numbers are reproduced by `literature_compare.py` (in this folder). It reads `results_fine.csv` only and does no new calculation. Run it with `python3 paper1/bo_offset/massfunc/literature_compare.py`. Statements marked "(from memory, unverified)" were not checked against a source.

## 0. Our convention (from THEORY.md)

* The nuclear reduced mass μ_N = M/2 becomes the tensor μ_N + 2A(R). Here A_par = ⟨∂_Rφ0|Q(h−e0)^{-1}Q|∂_Rφ0⟩ (Σg states only) and A_perp = Σ_n |⟨0|L_x|n⟩|²/(R²(e_n−e0)) (Πg states only), with unit electron mass.
* m_v = 4A_par and m_r = 4A_perp, entering as M → M + m. Equivalently, to leading order in 1/M, μ_v = (M + m_v)/2 and μ_r = (M + m_r)/2. The argument is r = μ_e R. At leading order in 1/M this equals R, so the scaling only matters at O(1/M²) in the energy.
* Separated-atom limit m_v = m_r = 1/2. The reduced mass of X+ + (X e) is M(M+1)/(2M+1) ≈ (M + 1/2)/2, so each nucleus carries half an electron mass. For H2 (two electrons, one per atom) the same counting gives 1.

## 1. Bunker and Moss 1977 [BM]

**What it computes.** It derives the effective vibration-rotation Hamiltonian for a ¹Σ diatomic to second order. The contact transformation is done without treating W_n(R) as constants. Nuclear-centre-of-mass coordinates are used, with electron mass m in H0 and mass polarization (1/2M_N)P² in H′. The result is Eq. (11),

    H_eff = (1/2μ)[P_R² + P_R β(R) P_R] + W̃(R) + B(1 + α(R)) J²,

with (Eq. 12) β(R) = −2μ Σ_n B_0n²/Δ_0n, where B_0n = (1/μ)⟨0|(P_R|n⟩), the sum runs over ¹Σg+ states and Δ_0n = W0 − W_n < 0. The rotational term is α(R) = B Σ_n |⟨0|L_x|n⟩|²/Δ_0n, summed over ¹Πg states. The reduced masses are μ_v(R) = μ/(1+β), μ_r(R) = μ/(1+α) (p. 422). It also gives a nonadiabatic potential ΔW = ΔW1 + ΔW2. For heteronuclear isotopologues there are extra μ/μ_a² terms (Section 4). Everything is formal. No function is computed. The only numbers are order-of-magnitude values for H2 (two electrons): α ≈ −1.6×10⁻⁵ (from the rotational g factor) and β "negative, of order up to about 10⁻³" (Unsöld estimate). The summary also says "−10⁻³ and −2×10⁻⁵".

**Convention mapping.** With P_R = −i∂_R (a.u.), B_0n² = −|⟨0|∂_R n⟩|²/μ², so

    β(R) = −(2/μ) A_par(R),   μ_v = μ/(1+β) ≈ μ + 2A_par   →   **m_v = −Mβ = 4A_par**   (μ = M/2).

This is exactly our Eq. (3) and our "μ_N + 2A". The sign also agrees: β < 0, so the mass increases. Their H0 has electron mass 1 rather than μ_e. Mass polarization sits in W_ad and ΔW rather than in the electronic Hamiltonian. These differences change m_v only at relative order 1/M.

For the rotation term, our derivation gives second-order perturbation theory in −2B(J_xL_x + J_yL_y) as α = 4B Σ_n |⟨0|L_x|n⟩|²/Δ_0n = −(2/μ)A_perp, and hence **m_r = −Mα = 4A_perp**. The printed [BM] formula lacks this factor 4. Two checks favour our version:
1. With the printed α, μ_r → μ + A_perp(∞)/2 = μ + 1/16 at large R, instead of the exact μ + 1/4. Kutzelnigg's B(∞) = −1/2 (Section 3 below) gives μ + 1/4.
2. The printed α(H2) ≈ −1.6×10⁻⁵ corresponds to m_r = 0.029. Four times it gives m_r = 0.1175, which equals 1 − g_J(H2) for g_J(H2) ≈ 0.883 nuclear magnetons (g value from memory, unverified). In our second-order theory g_J = 1 − 4A_perp = 1 − m_r for a homonuclear molecule (our derivation, not in the papers).

The factor 4 may be a misprint or an unstated normalization of the Π states. Treat it as unresolved, but it does not affect β.

**Numerical comparison.** No H2+ values are given. For H2 (two electrons) we can compare only orders of magnitude. Our H2+ values at R_e are β = −m_v/M = −1.9×10⁻⁴ and α = −m_r/M = −4.1×10⁻⁵. These are the same orders as their H2 estimates, which is all one can say.

**Bearing on our claims.** [BM] establishes the formal object we compute: R-dependent μ_v(R), μ_r(R) from a second-order sum over Σg and Πg states. Our m_v, m_r are their β, α (up to −M and the factor-4 question for α), evaluated exactly for H2+. The idea and the formulas are therefore not new. What is new relative to [BM] is the evaluation.

## 2. Bunker, McLarnon and Moss 1977 [BMM]

**What it computes.** It fits [BM]'s Eq. (5) with **constant** μ_v, μ_r (R dependence ignored) to D2 levels with v ≤ 20 and J ≤ 5, using the Kołos–Wolniewicz W(R). μ_r − μ = 0.000014 ± 0.000034 u is not significant, so μ_r = μ is fixed. The fit gives μ_v = μ + 0.0002538(12) u (rms 0.096 cm⁻¹). Through the isotope scaling (μβ isotope invariant, μβ = −0.0002537 u) they predict H2: μ_v = 0.5038920 u, against μ = 0.5036382 u and μ_atomic = 0.5039125 u. They also give Unsöld (average-energy) estimates of ⟨0|∂²/∂R²|0⟩ = −0.47 Å⁻², ⟨L²⟩ = 0.37 and α(D2) = −4.5×10⁻⁵. They state the bracketing μ < μ_v < μ_atomic as the physical picture (non-adiabatic lag).

**Convention mapping and numbers** (two electrons, atomic limit m = 1):

| quantity | [BMM] | in our units (m_e, M → M + m) |
|---|---|---|
| fitted constant μ_v − μ (D2, v ≤ 20) | 0.0002538 u | m_v,eff = 0.925 (92.5% of the atomic value 1.000) |
| fitted μ_r − μ | 0.000014 ± 0.000034 u | m_r = 0.05 ± 0.12 |
| Unsöld α(D2) | −4.5×10⁻⁵ | m_r = 0.165 (printed formula) or 0.66 (with factor 4) |

**What can be compared.** They treat H2 and D2 only, with two electrons, so there is no direct number for H2+. One qualitative parallel holds. A single constant vibrational mass fitted to many levels comes out close to the atomic value (92.5%), even though theory gives a much smaller value near R_e (Kutzelnigg: 65% for H2 at R_e). This is the same pattern we find for H2+. Our m_v(R_e) is 0.347 (69%), but the threshold states give ⟨m_v⟩_w = 0.48–0.50 (96–99%), because highly excited levels sample the large-R region and the maximum of m_v(R). [BMM] state that the systematic residuals "could no doubt be rectified by introducing an R-dependence into μ_v". Our functions supply that R dependence for H2+. The insensitivity of their fit to μ_r is consistent with a small m_r near R_e (ours: 0.075, 15% of the limit).

## 3. Kutzelnigg 2007 [K]

**What it computes.** He builds a non-adiabatic ansatz that is exact in every dissociation channel (Sections 6.4–6.5), Ψ = ψ_gχ + (1/(M+1)) ω·∇χ, and truncates the energy at O(M⁻²). The result is a single potential with geometry-dependent reduced masses (Eq. 120):

    1/μ_vib = (1/M_R)(1 + A(R)/M),   1/μ_rot = (1/M_R)(1 + B(R)/M),   M_R = M/2,

with A = −4⟨ω̃_z|∂ψ_g/∂Z⟩ + ⟨ω̃_z|H_cn − E_cn|ω̃_z⟩ (Eq. 121) and B the same with x (Eq. 122). This is a Hylleraas functional evaluated with a fixed trial ω (the asymptotically exact one) rather than its minimizer. In Section 6.9 he writes the stationarity condition (H_cn − E_cn)τ = −∇ψ (Eq. 151), which is the Sternheimer equation, but he does not solve it. Conclusion (a) lists "formulate the theory in terms of exact wave functions" as future work. Numbers are LCAO for H2+ and Heitler–London for H2.

**Convention mapping.** μ_vib ≈ (M − A)/2, so **m_v = −A and m_r = −B**. At the stationary point, −A = 4⟨∂_Rφ0|(H−E)^{-1}|∂_Rφ0⟩ = 4A_par, which is identical to our definition (similarly −B = 4A_perp). His asymptote A(∞) = B(∞) = −1/2 (Eq. 135), with 1/M_R(1 + A/M) = 1/M_ρ (Eq. 136), equals our m(∞) = 1/2. He uses the clamped-nuclei electronic Hamiltonian with electron mass 1, so **a clamped-nuclei convention does not double the asymptote**. The statement in THEORY.md §5 (last bullet) and in the draft Section 6.5 ("a clamped-nuclei convention doubles the asymptote") appears to be wrong. A value of 1 arises only if one counts the correction per total nuclear mass 2M, or for H2 with two electrons. "Participation" in [K] = m/(1/2).

**Numerical comparison for H2+** (values read from the text of [K]; figure-read values marked ≈):

| R (bohr) | [K] model | −A_K | our m_v | −B_K | our m_r |
|---|---|---|---|---|---|
| 2.0 | simple LCAO, α = 1 | 0.20 | 0.3467 | 0.06 | 0.0756 |
| 2.5 | simple LCAO, α = 1 | 0.29 | 0.3804 | 0.11 | 0.1074 |
| 2.0 (R_e) | Finkelstein–Horowitz, optimized α | **0.31** | **0.3467** | **0.14** | **0.0756** |
| ≈3 | simple LCAO: A back to −0.5 | 0.5 | 0.4254 (ours reaches 0.5 at R = 3.68) | – | 0.1436 |
| ≈5 | simple LCAO: minimum A ≈ −0.6 | 0.6 | 0.6317 (ours: maximum 0.6606 at R = 5.905) | – | 0.3114 |
| → ∞ | both | 0.5 | 0.5 (as R⁻³) | 0.5 | 0.5 (as R⁻³) |
| → 0 | both | 0 | 0 (∝ R^1.87 at R = 0.02) | 0 | 0 (∝ R²) |

Agreement and disagreement:
* **Asymptotes and topology agree exactly.** Both A and B go to −1/2 at large R and to 0 at R = 0. A(R) overshoots the atomic value at intermediate R, which is the "hump". He writes that there is a region where "the participation of the electron is slightly enhanced" (p. 2643). B(R) is monotonic, and he attributes this to ω_z being a bonding pσ combination while ω_x is an antibonding π combination. Our m_v maximum (0.661 at 5.9 bohr) matches his simple-LCAO minimum (≈ −0.6 near 5 bohr) in size and position. His Finkelstein–Horowitz curve (Fig. 6) has its minimum near ρ = αR ≈ 4, roughly R ≈ 3.5 (figure-read, approximate).
* **At R_e his best (FH) values are off by −11% (vibration) and +86% (rotation).** He gives 0.31 against our 0.3466 (60% against 69% participation) and 0.14 against our 0.0754 (about 30% against 15%). The rotational mass is therefore much less "nuclear" in his LCAO estimate than in the exact second-order result. His slogan "60% in vibration, 30% in rotation" becomes about 70% and 15% for H2+. Our values are the exact stationary points of the Hylleraas functional, which he uses with a fixed trial function. They also agree, untuned, with the values inferred from the exact three-body energies (m_v from c6: 0.3465 ± 0.0003; m_r from the J spacings: 0.0754). So the discrepancy is in [K]'s model wave functions. For H2 he gives 65% and 35% (A ≈ −0.65, B ≈ −0.37, with A ≈ A1/2). If g_J(H2) ≈ 0.883 (memory, unverified), the exact second-order m_r(H2) ≈ 1 − g_J ≈ 0.12, i.e. about 12% participation. That suggests the same LCAO/HL overestimate of the rotational share.
* **Introduction of [K], prior H2+ work cited there.** R. E. Moss, Mol. Phys. **89**, 195 (1996) interpreted non-adiabatic H2+ calculations "a posteriori in terms of effective moving masses". The message was atomic masses for vibration and nuclear masses for rotation as good compromises. The same recipe was used for H3+ (Röhse, Kutzelnigg, Jaquet, Klopper, JCP **101**, 2231 (1994), and follow-ups) and for HF (Müller et al., Theor. Chem. Acc. **100**, 85 (1998)). Herman and Asgharian, J. Mol. Spectrosc. **19**, 305 (1966) and Fisk and Kirtman, JCP **41**, 3516 (1964) are earlier effective-Hamiltonian work (all bibliographic data as printed in [K] and [BM]). Moss 1996 is the closest prior numerical H2+ work we know of from these papers and has **not** been read.

## 4. Novelty assessment (honest)

* **Not new:**
  * The concept and the second-order formulas for R-dependent vibrational and rotational masses ([BM] 1977; earlier Herman–Asgharian 1966, Fisk–Kirtman 1964).
  * The identification of their large-R limit with the atomic (nucleus + m_e/2 per nucleus for H2+) reduced mass ([K] Eqs. 135–136; [BMM] bracketing μ < μ_v < μ_atomic).
  * The qualitative shape for H2+: a vibrational hump above the atomic value, a monotonic and small rotational mass, both vanishing at R → 0 ([K] Figs. 4, 6, 7, LCAO level).
  * The exact defining equation (our Sternheimer equation is [K] Eq. 151, stated there but not solved).
* **Probably not new either (from memory, unverified; check before submission):** accurate R-dependent nonadiabatic mass functions for H2+ and HD+ have, we believe, been computed after 2007. Candidates:
  * R. Jaquet and W. Kutzelnigg, Chem. Phys. **346**, 69 (2008), "Non-adiabatic theory in terms of a single potential energy surface. The vibration–rotation levels of H2+ and D2+";
  * D. W. Schwenke, J. Chem. Phys. **114**, 1693 (2001), on an effective Hamiltonian for H2+ and HD+;
  * K. Pachucki and J. Komasa, J. Chem. Phys. **129**, 034102 (2008) and **130**, 164113 (2009), H2;
  * E. Mátyus and S. Teufel, J. Chem. Phys. **151**, 014113 (2019), mass-correction functions;
  * Moss 1996 (cited in [K]).

  All citation details here are from memory except Moss 1996. If any of these tabulates exact m_v(R), m_r(R) for H2+, our functions are a recomputation. Theirs should then be compared numerically, and ours cited as an independent check.
* **What is ours:**
  1. A converged evaluation (3×10⁻¹¹) from the exact BO wave function over 0.02–60 bohr, including the R⁻³ approach to 1/2 (coefficients 18.19 and −9.03) and the small-R power laws. We have not seen these in the three papers.
  2. Their use, with no fitted parameter, to explain the near-constant BO offsets of the threshold critical masses and the κ⁶ coefficient. That application is not in any of the three papers.
  3. The untuned agreement with m_v(R_e) and m_r(R_e) extracted from exact three-body energies.
* **"First-principles".** The label is defensible only in the sense the draft already gives ("from the exact BO solutions, no fitted parameter, second order"). It should not suggest a new theory. Better: "second-order nonadiabatic (Bunker–Moss) mass functions evaluated from the exact BO wave function".
* **"Computed for the whole family".** In the μ_e-scaled convention the unit-mass functions m(r) do not depend on M, so one function serves every member. This is a scaling property shared by any H2+ calculation, not a separate computation, and it holds only to the order of the theory (the neglected terms are O(1/M) in λ_c). Say "one universal function serves the whole family" rather than "computed for the whole family".

## 5. Suggested replacement wording (plain prose)

**Draft Section 6.2, paragraph "First-principles mass functions" (replace the heading and first sentences):**

> **Non-adiabatic mass functions (v1.6, `paper1/bo_offset/massfunc/`).** The effective vibration-rotation Hamiltonian of Bunker and Moss [109] contains R-dependent vibrational and rotational reduced masses given by second-order sums over Σg and Πg electronic states. Kutzelnigg [109] estimated them for H₂⁺ with LCAO wave functions and wrote down, but did not solve, the exact defining equation. We evaluated these functions from the exact BO wave function by solving that equation (a Sternheimer equation, prolate-spheroidal basis, 2597 values of R from 0.02 to 60 bohr, basis convergence 3×10⁻¹¹). In the μ_e-scaled convention a single pair of functions m_v(R), m_r(R) serves every member of the family. m_v(R) rises from 0.347 at R_e to a maximum of 0.661 at R = 5.9 bohr and returns to 1/2 as R⁻³. m_r(R) rises monotonically from 0 to 1/2. The limit 1/2 is the reduced mass of X⁺ and the atom. These shapes agree with Kutzelnigg's LCAO estimates, which also show the vibrational maximum above 1/2 and a monotonic rotational mass. At R_e our values are 12% larger (vibration) and 46% smaller (rotation) than his best estimates (0.31 and 0.14). At R_e our values reproduce, with no adjustment, the values inferred from the exact energies, m_v(R_e) = 0.34655 (from c₆, Section 6.5: 0.3465 ± 0.0003) and m_r(R_e) = 0.07542 (from the J spacings: 0.0754).

and replace the closing clause of that paragraph ("The next-order term, ... and the comparison with published effective vibrational and rotational mass functions (literature to be checked) are outstanding, so "first-principles" here means ...") by:

> The next-order term and the sensitivity to the mass convention (Section 6.5) are outstanding. Accurate mass functions for H₂⁺ may exist in later work that we have not yet checked. The new element here is not the functions themselves but their use, with no fitted parameter, to account for the threshold offsets and the κ⁶ coefficient.

**Abstract clause** (currently "...are accounted for, with no fitted parameter and at the 1–2% level of the offset for v=1–3, by non-adiabatic vibrational and rotational mass functions computed to second order in perturbation theory (not yet compared with published mass functions)"):

> ...are accounted for, with no fitted parameter and at the 1–2% level of the offset for v = 1–3, by the second-order (Bunker–Moss) non-adiabatic vibrational and rotational mass functions, which we evaluate from the exact Born–Oppenheimer wave function.

If the shorter wording "follow, with no fitted parameter, from first-principles non-adiabatic mass functions" is used anywhere, replace it with "follow, with no fitted parameter, from second-order non-adiabatic mass functions evaluated from the exact Born–Oppenheimer wave function".

**Draft Section 6.5, last sentence** ("We have not compared m_v(R) with published vibrational-mass functions for H₂⁺ (not checked; conventions differ, e.g. a clamped-nuclei convention doubles the asymptote).") should become:

> In Kutzelnigg's convention [109] (clamped-nuclei electronic Hamiltonian, 1/μ_vib = (2/M)(1 + A/M)), m_v = −A, with the same asymptote 1/2; his LCAO estimate at R_e is 0.31.

**THEORY.md §5, last bullet** should be corrected in the same way.

**Reference [109]** is already in the draft with [BM], Herman–Asgharian and [K]. Suggested additions, with data as printed in the PDFs read:
* P. R. Bunker, C. J. McLarnon, R. E. Moss, Mol. Phys. **33**, 425 (1977), "Application of the effective vibration-rotation hamiltonian to H₂ and D₂", doi:10.1080/00268977700100361.
* R. E. Moss, Mol. Phys. **89**, 195 (1996) (as cited in [K]; title not seen; to be read).
* G. A. Fisk, B. Kirtman, J. Chem. Phys. **41**, 3516 (1964); D. M. Bishop, S.-K. Shih, J. Chem. Phys. **64**, 162 (1976) (as cited in [BM]).

## 6. Open items

* Read Moss 1996 and check the memory-only candidates in Section 4 (especially Jaquet and Kutzelnigg 2008 and Schwenke 2001) for tabulated H₂⁺ mass functions. Compare numerically if they exist.
* Resolve the factor 4 in [BM]'s α (it affects only the rotational comparison).
* Optional cross-check that needs no new heavy calculation: g_J(H₂⁺) ≈ 1 − ⟨m_r⟩ ≈ 0.92–0.93 (our relation) against a published H₂⁺ rotational g factor, if one exists.

## 7. Moss 1996 [M96]

[M96] R. E. Moss, Mol. Phys. **89**, 195–210 (1996), "On the adiabatic and non-adiabatic corrections in the ground electronic state of the hydrogen molecular cation", doi:10.1080/002689796174083. Read in full ( Table 5 was read from a rendered page image because the text layer is rotated and garbled).

**What it computes.** [M96] computes **no mass functions**. It takes the exact non-adiabatic dissociation energies of all H₂⁺ levels (Moss, Mol. Phys. **80**, 1541 (1993): 462 levels, v ≤ 19, N ≤ 40, including quasibound ones; D₂⁺: 164 levels with N ≤ 5) and subtracts the "partitioned adiabatic" energies. This gives the non-adiabatic correction for every level (Table 4, printed to 0.001 cm⁻¹; maximum 1.3 cm⁻¹ near v = 10–11). It then fits the Bunker–Moss effective equation (his Eq. 11) to these numbers with **constant** μ_v, μ_r and a potential scale factor (1 + γ) (Tables 5, 6). He states that attempts to make the parameters R dependent failed because the extra parameters were highly correlated (p. 206). His closing suggestion (p. 208) is that "a suitable R dependence for the reduced masses" might remove the need for potential scaling. So R-dependent μ_v(R), μ_r(R) are named in [M96] as open, not computed. Other content: the adiabatic correction for every level (Table 1), a scaled-BO-potential model for it (Table 3), a critique of the Röhse et al. reduced-electron-mass recipe, and the isotope relations μ_v = μ_nuc/(1 + b/μ_nuc), μ_r = μ_nuc/(1 + a/μ_nuc) (Eqs. 13–14) tested by predicting D₂⁺ from H₂⁺.

**Convention mapping.** [M96] quotes masses in u with μ_nuc = 0.503638235 u and μ_at = 0.503775380 u (half an electron mass added per nucleus). His "participation" p = (μ − μ_nuc)/(μ_at − μ_nuc) (Table 8) equals **2m** in our units (our atomic limit is m = 1/2). His isotope constants map as b = −m_v·(m_e in u)/2 and a = −m_r·(m_e in u)/2, i.e. p_v = −b/(μ_at − μ_nuc), p_r = −a/(μ_at − μ_nuc). His reference for "non-adiabatic correction" is the exact adiabatic potential with nuclear reduced mass, which is the same reference as our BO+ad with M → M + m. Our m(r) uses r = μ_e R; at M = 1836 this differs from R by 3×10⁻⁴ relative and is ignored below.

**His fitted constants in our units** (from `literature_compare.py`):

| [M96] fit | p_v | m_v,eff | p_r | m_r,eff |
|---|---|---|---|---|
| Table 5 model 6 (μ_v, μ_r; γ = 0; 462 levels) | 0.828 | 0.414 | 0.343 | 0.172 |
| Table 5 model 9 (μ_v, μ_r, γ; 462 levels) = Table 8 | 0.951(8) | 0.475 | 0.416(5) | 0.208 |
| Table 6 H₂⁺ N ≤ 5 (all three) | 0.891(7) | 0.445 | 0.28(10) | 0.14 |
| Table 7 isotope constants b, a | 0.890 | 0.445 | 0.28 | 0.14 |
| ours at R_e | 0.693 | 0.3466 | 0.151 | 0.0754 |
| ours, ⟨m_v⟩ weighted by \|χ′\|² for v = 0 … 19 (N = 0) | 0.70 … 0.97 | 0.351 … 0.486 | | |
| ours, constant-mass least-squares fit to our own predicted corrections (423 bound levels, N ≤ 40, γ = 0; analogue of model 6) | 0.840 | 0.420 | 0.307 | 0.154 |

His constants are level-averaged effective values that depend on which parameters are freed (p_v ranges from 0.77 to 1.03 across models 3–9). They are not values of a function at any R. The analogue of his model 6 fitted to our corrections gives 0.840 and 0.307 against his 0.828 and 0.343 (rms 0.109 against his 0.135 cm⁻¹; he also includes 39–40 quasibound levels). His conclusion "atomic masses for vibration, nuclear masses for rotation" is the level-averaged shadow of m_v rising from 0.35 to 0.49 across the vibrational ladder while m_r stays at 0.08–0.19 for N ≤ 35.

**Direct numerical test: his Table 4 from our functions, with no fit.** To first order in m(R)/M, with BO vibrational functions χ_vN,

    ΔD(v,N) = (1/M²) [ ∫ m_v(R) χ′² dR + N(N+1) ∫ m_r(R) χ²/R² dR ],   M = 1836.15267 m_e.

`literature_compare.py` evaluates this by a 1D finite-difference solve of E0_exact(R) + 1/R from `results_fine.csv` (h = 0.002 bohr, box 0.4–60 bohr; going from h = 0.004, box 40 changes the corrections by about 0.001 cm⁻¹ or less, except v = 19, 0.7 cm⁻¹ below threshold: 0.012 → 0.009). The BO D₀ comes out as 21375.93 cm⁻¹ with 20 bound N = 0 levels, as expected. Result against all 462 printed entries of Table 4:

* 423 entries refer to levels bound in our BO potential; 39 quasibound ones are skipped.
* Two entries differ grossly and are evidently misprints in [M96]: v = 0, N = 2 printed 0.177 (ours 0.116; neighbours 0.111 and 0.123) and v = 10, N = 4 printed 1.037 (ours 1.308; neighbours 1.309 and 1.305).
* For the other 421: **max |ours − Moss| = 0.0036 cm⁻¹, rms 0.0010 cm⁻¹**, mean +0.0007, median 0.0007 cm⁻¹. 95.5% agree to ≤ 0.002 cm⁻¹, i.e. within the printed rounding (0.0005) plus about 0.001. The N = 0 column agrees to ≤ 0.002 cm⁻¹ for every v from 0 to 19, including the turnover at v = 11.
* Rotational part, checked separately: at v = 0, N = 35 our correction is 0.062 (vibration) + 1.211 (rotation) = 1.273 against his 1.271 cm⁻¹. A rotational function four times smaller or larger would miss by about 0.9 or 3.6 cm⁻¹.

So the second-order mass functions alone reproduce the exact non-adiabatic corrections of the whole H₂⁺ vibration-rotation manifold to his printed precision. The remaining O(M⁻²) pieces (the Bunker–Moss potential ΔW, the R → μ_e R argument scaling, adiabatic corrections to χ) together shift dissociation energies by at most a few 10⁻³ cm⁻¹. Caveats: this is a first-order expectation value with BO χ, not a solution of the effective equation, and Table 4 was parsed from the PDF text layer (spot-checked against the page, not independently retyped).

**Factor-4 question (Section 1).** [M96] does not restate the α formula, so it does not settle the [BM] printing. Numerically, however, the question is settled for our purposes: m_r = 4A_perp reproduces the N dependence of his exact corrections to 0.002 cm⁻¹, whereas the printed [BM] normalization (m_r = A_perp) would give a rotational contribution four times too small.

**Novelty verdict (honest).** Relative to [M96] our work is an **extension, not a reproduction**. Moss had the exact level-by-level non-adiabatic corrections, extracted only constant effective masses from them, explicitly failed to determine an R dependence, and named R-dependent reduced masses as future work. We evaluate the R-dependent Bunker–Moss functions directly from the exact BO wave function and show that, with no fitted parameter, they reproduce his entire Table 4. This also shows that the R dependence he could not fit is fully determined by second-order theory. What is still not ours: the formalism ([BM], Herman–Asgharian, Fisk–Kirtman), the qualitative shapes ([K]), the exact non-adiabatic energies (Moss 1993 and later). The remaining risk is post-1996 work not yet read, especially Jaquet and Kutzelnigg 2008 and Schwenke 2001 (Section 4; from memory, unverified). Either may tabulate accurate H₂⁺ mass functions. [K] (2007) does not cite any such tabulation and still calls an exact-wave-function treatment future work, which weakly suggests none existed before 2007.

**Suggested wording.**

Replace "Accurate mass functions for H$_2^+$ may exist in later work that we have not yet checked (for example Moss [112]). The new element here is not the functions themselves but their accurate evaluation and their use..." by:

> Moss [112] fitted constant effective vibrational and rotational masses to the exact non-adiabatic level energies of H$_2^+$ and reported that an R dependence could not be determined from the energies. Our functions supply that R dependence from theory alone. Inserted to first order with Born-Oppenheimer vibrational functions, they reproduce all 421 non-misprinted entries of his table of non-adiabatic corrections (v ≤ 19, N ≤ 40, corrections up to 1.3 cm$^{-1}$) to within 0.004 cm$^{-1}$ (rms 0.001 cm$^{-1}$), with no fitted parameter. His preferred constants (vibrational participation 0.83 to 0.95, rotational 0.3 to 0.4 of the atomic value) are level averages of our functions. The formalism is that of Bunker and Moss [109], and we have not yet checked later work (for example Jaquet and Kutzelnigg 2008) for tabulated functions. The new elements here are the accurate evaluation of the functions over the whole range of R, this parameter-free test against the exact H$_2^+$ manifold, and their use...

Replace the ref [112] note "(verified via Crossref; not yet read)" by:

> R. E. Moss, Mol. Phys. 89, 195 (1996), doi:10.1080/002689796174083 (read in full; constant effective masses fitted to exact non-adiabatic H$_2^+$ and D$_2^+$ levels; no R-dependent functions).

**Items in Sections 3, 4 and 6 now superseded.** "Moss 1996 ... has not been read" (Section 3) and "Read Moss 1996" (Section 6) are done. The factor-4 item in Section 6 is resolved numerically (above). Still open: check Jaquet–Kutzelnigg 2008 and Schwenke 2001.

## 8. Schwenke 2001 [S01]

[S01] D. W. Schwenke, "A first principle effective Hamiltonian for including nonadiabatic effects for H2+ and HD+", J. Chem. Phys. **114**, 1693–1699 (2001), doi:10.1063/1.1334897. Read in full ( Fig. 1 read from a 300 dpi render, Table I from the text layer and checked on the page).

**What it computes.** Exactly the object we compute, for H₂⁺ and HD⁺, and it is the first paper to do so (his words: to his knowledge the explicit Bunker–Moss expressions "have not been used to generate the effective Hamiltonian"). He keeps the full second-order Bunker–Moss operator without removing the first-derivative term (Eq. 15):

    H_na = (∂_R)† C_vv ∂_R + [(∂_R)† C_v + C_v ∂_R] + Σ_{a=x,y} N_a† C_rr N_a + C,

with C_vv = Σ_n b_v^n b_v^n/Δ_0n over Σ states (Eq. 16), C_rr = Σ_n b_x^n b_x^n/Δ_0n over Π states (Eq. 18), and the cross term C_v and potential term C (Eqs. 17, 19). These contain the mass-polarization and ⟨n|L_a|0⟩-type pieces in b^n (Eq. 11).

* **Method.** Electronic states from MOLPRO with nuclear-centred Gaussian (aug-)cc-pVnZ basis sets, all excited states of the basis included in the sums (a discrete sum over a finite basis, not a Sternheimer solve). Derivative couplings by symmetric finite differences at ΔR = 0.05 and 0.025 bohr with Richardson extrapolation. Augmented basis sets were needed for stability at R > 4 bohr. No convergence estimate is given for the functions themselves. The calculated R grid is not stated. The functions were interpolated (R⁶-weighted Lagrange) and extrapolated by straight lines below and decaying exponentials above the data. Fig. 1 shows R = 1–20 bohr.
* **Output for the functions.** Figures only (Fig. 1a–e: adiabatic correction, C_vv, C_v, C, C_rr R² for H₂⁺ and HD⁺). The numbers are in EPAPS E-JCPSA6-114-005105 (one ASCII file of 15 kB: "the ab initio calculations of the nonadiabatic correction for H2+, HD+"). **Not seen.** It almost certainly tabulates the functions on his R grid.
* **Level test.** He solves the effective equation variationally (Laguerre basis and a finite-difference/trapezoid expansion, R = 0.2–150 bohr) on the Bishop–Wetmore/Peek BO potential plus the Bishop–Wetmore adiabatic correction, with Moss's masses, for every bound and long-lived quasibound level with N = 0, 10, 20, 30 (Table I). Against Moss's exact non-adiabatic energies, the dissociation-energy errors are 10.6 (BO), 1.31 (adiabatic) and **0.00112 cm⁻¹ max (non-adiabatic)**, with rms 6.9, 1.00 and **0.00062 cm⁻¹**. ΔG₀ is off by −0.00007 cm⁻¹, and most transition frequencies agree to about 10⁻⁴ cm⁻¹. HD⁺ is worse (rms 0.009 cm⁻¹, v ≤ 19), which he attributes to his own adiabatic correction. In HD⁺ the functions blow up at 9–15 bohr, where the electron localizes on D (g/u mixing). This breaks the second-order treatment for v ≥ 20. H₂⁺ has no such problem.
* **Comment on Moss.** He notes that Moss's constant μ_v, μ_r correspond to constant C_vv and C_rr R², and that his functions agree with Moss's values near R ≈ 3 bohr (our m(3) = 0.425, 0.144 against Moss model 6: 0.414, 0.172).

**Convention mapping.** G_v = −1/(2μ_R), μ_R = M/2 and b_v^n = 2G_v⟨n|∂_R 0⟩ give C_vv = −A_par/μ_R² = **−m_v/M²**. With G_r = 1/(2μ_R R²) and b_x^n = −2G_r⟨n|L_x|0⟩, C_rr = −A_perp/(μ_R² R²) = **−m_r/(M² R²)**. The plotted quantities are therefore C_vv×10⁷ = −0.33715 m_v and C_rr R²×10⁹ = −337.15 m_r (M = 1836.152701). The asymptote m = 1/2 gives −1.483 and −148.3. The figure shows about −1.49 and −148, which confirms the mapping, the sign and the factor 4 in m_r (Section 1).

**Numerical comparison** (`literature_compare.py`, digitized from Fig. 1b and 1e; the digitization uncertainty is about 1% of full scale, i.e. ±0.005 in m_v and ±0.003 in m_r):

| R (bohr) | m_v [S01] | m_v ours | m_r [S01] | m_r ours |
|---|---|---|---|---|
| 1.5 | 0.323 | 0.3233 | 0.047 | 0.0484 |
| 2.0 | 0.347 | 0.3467 | 0.074 | 0.0756 |
| 3.0 | 0.424 | 0.4254 | 0.141 | 0.1436 |
| 4.0 | 0.536 | 0.5372 | 0.223 | 0.2262 |
| 6.0 | 0.658 | 0.6603 | 0.380 | 0.3828 |
| 8.0 | 0.590 | 0.5917 | 0.460 | 0.4616 |
| 10.0 | 0.536 | 0.5357 | 0.486 | 0.4865 |
| 19.5 | 0.503 | 0.5025 | 0.499 | 0.4988 |

Over 17 points from 1.5 to 19.5 bohr: m_v differs by mean −0.0005 and max 0.0026, and m_r by mean −0.0016 and max 0.0035. **His functions and ours agree to the precision of the figure.** That includes the maximum of m_v near 6 bohr and the slow approach to 1/2. The small systematic −0.002 in m_r is within the digitization error (the digitized C_rr asymptote reads 0.4987). A sharper comparison needs the EPAPS tables.

**His level test against ours.** Both reproduce Moss's exact H₂⁺ non-adiabatic energies to about 10⁻³ cm⁻¹. His test is stronger in method: he solves the full effective equation, including C_v and C, with an independent BO+adiabatic potential, and gets max 0.0011 and rms 0.0006 cm⁻¹ on total dissociation energies for N = 0, 10, 20, 30. Our test (Section 7) is a first-order expectation value of the mass terms alone on Moss's Table 4 corrections, with max 0.0036 and rms 0.0010 cm⁻¹ over N ≤ 40. Ours covers more levels (all N ≤ 40) and shows that the mass functions alone already carry the corrections at this precision. His Fig. 1 shows C and C_v that are not individually negligible: C(R) − C(∞) ≈ +0.5×10⁻⁷ E_h ≈ 0.01 cm⁻¹ near R_e. In his partition (Bunker–Moss without the first-derivative elimination, mass polarization in H′) they must largely cancel against each other or against the difference between his and our partition of adiabatic and non-adiabatic terms. We have not resolved this (open item). The statement "the remaining O(M⁻²) pieces shift energies by at most a few 10⁻³ cm⁻¹" (Section 7) holds in our partition as an empirical result, not as a statement about each of his terms.

**Bearing on our claims.** [S01] removes two of the three "new elements" claimed in the draft:
1. "Accurate evaluation of the functions over the whole range of R" is not new. He computed accurate ab initio C_vv(R), C_rr(R) (plus C_v, C) for H₂⁺ and HD⁺ in 2001 and published them as figures and EPAPS tables. Our evaluation is more precise and better documented: an exact BO wave function rather than a Gaussian basis, a direct Sternheimer solution with stated 3×10⁻¹¹ basis convergence, 0.02–60 bohr, analytic R⁻³ tail coefficients and small-R power laws. But the functions were known to graphical accuracy, and, as his 6×10⁻⁴ cm⁻¹ level agreement shows, to better than needed for spectroscopy.
2. "A parameter-free test against the exact H₂⁺ manifold" is not new. He did it in 2001, with a more complete Hamiltonian, and found the same 10⁻³ cm⁻¹ agreement with Moss. Our Section 7 statement "Relative to [M96] our work is an extension ... the R dependence he could not fit is fully determined by second-order theory" is true, but [S01] showed it first (and [S01] cites [M96]).

## 9. Jaquet and Kutzelnigg 2008 [JK08]

[JK08] R. Jaquet, W. Kutzelnigg, "Non-adiabatic theory in terms of a single potential energy surface. The vibration–rotation levels of H2+ and D2+", Chem. Phys. **346**, 69–76 (2008), doi:10.1016/j.chemphys.2008.02.068. Read in full ( Tables 4–6 checked on the rendered pages, because the text layer drops all minus signs).

**What it computes.** It implements Kutzelnigg's 2007 single-surface theory ([K], Section 3). The functions A(R) and B(R) of Eqs. (42)–(43) (same as [K] Eqs. 121–122) are evaluated with a **fixed** trial ω̃ built from model wave functions: simple LCAO (α = 1), Finkelstein–Horowitz (FH, LCAO with optimized exponent), and for B also FH plus a p polarization function (B(pol)). They also try two semi-empirical rescalings of B (M4, M5). They state that the results with FH "were so satisfactory that we did not care to go beyond this approximation", and that the best strategy for A and B "has probably not yet been found". The functions are closed-form LCAO expressions, **tabulated in Table 6** (R = 0.2–20 and 50 bohr, 86 points: A, B, B(pol) for H₂⁺, all printed negative). Nuclear motion uses the Kołos BO+adiabatic potential and a DVR grid (about 350 points, 0.1–50 bohr). Level tests: J = 0, 1 of H₂⁺ and D₂⁺ against Hilico et al. (2000) exact energies (Tables 1–3), and J = 0, 10, 20, 30 of H₂⁺ against Moss's non-adiabatic corrections [M96] (Tables 4–5).

They describe this as "the first implementation of the new scheme" and say related ideas "have never led to a computational scheme for the calculation of energy levels with non-adiabatic accuracy, except in an entirely empirical context". They **do not cite [S01]**, although [S01] had done the ab initio Bunker–Moss calculation for H₂⁺ seven years earlier.

**Convention mapping.** As in [K]: 1/μ_vib = (2/M)(1 + A/M), so **m_v = −A and m_r = −B**. The asymptote is −0.5 (Table 6 at 50 bohr: −0.500000 for all three columns). Their Table 6 is model-dependent and is not the stationary (exact second-order) value.

**Numerical comparison of the functions** (Table 6 against ours, `literature_compare.py`):

| R | −A (FH) | m_v ours | −B (FH) | −B(pol) | m_r ours |
|---|---|---|---|---|---|
| 1.0 | 0.103 | 0.303 | 0.063 | 0.036 | 0.026 |
| 2.0 | 0.308 | 0.347 | 0.146 | 0.098 | 0.076 |
| 3.0 | 0.477 | 0.425 | 0.227 | 0.173 | 0.144 |
| 4.0 | 0.574 | 0.537 | 0.306 | 0.259 | 0.226 |
| 5.0 | 0.600 | 0.632 | 0.377 | 0.344 | 0.311 |
| 6.0 | 0.586 | 0.660 | 0.430 | 0.413 | 0.383 |
| 8.0 | 0.537 | 0.592 | 0.482 | 0.479 | 0.462 |
| 10.0 | 0.512 | 0.536 | 0.495 | 0.495 | 0.487 |

* FH vibrational maximum: 0.6003 at R = 5.0, against ours 0.6606 at 5.905. Between R = 1 and 6 the largest difference is 0.20 in m_v (at R = 1, where FH quenches far too strongly), 0.084 in m_r(FH) and 0.033 in m_r(pol).
* At R = 2: m_v(FH) is 11% low. m_r(FH) is 1.93 times ours and m_r(pol) 1.30 times ours. The polarization function moves B toward the exact value. Their remark that B, not A, is the critical function, and that it needs more quenching at intermediate R (4–6 bohr), is exactly the direction of the exact result.

**Their level tests, reproduced.** Using their Table 6 functions in our first-order formula (BO χ from our potential, Section 7) reproduces their printed deviations from Moss for M1 (FH A, B; 52 levels bound in our BO potential, J = 0, 10, 20, 30; J = 20, v = 11 and J = 30, v = 5 are quasibound and skipped) with rms 0.0027 cm⁻¹. The largest mismatch, 0.010 cm⁻¹, is at J = 0, v = 19, 0.7 cm⁻¹ below threshold. For M3 (FH A, B(pol); 32 levels) the rms is 0.0022 cm⁻¹. The remaining near-constant offset of about 0.002–0.003 cm⁻¹ reflects their DVR, the Kołos potential and the fact that they solve rather than take first order. This independently confirms the mapping m = −A, m = −B and our first-order procedure. Their errors (M1: −0.006 to −0.035 cm⁻¹ at J = 0, up to +0.52 cm⁻¹ at J = 30; M3: up to 0.19 cm⁻¹ at J = 30) therefore come entirely from the model functions. On the same 52 levels our functions give rms 0.0010 and max 0.0025 cm⁻¹, i.e. 10 to 200 times smaller errors. Against the absolute Hilico energies (J = 0, 1) their best H₂⁺ errors are below 0.06 cm⁻¹ and their D₂⁺ errors below 0.03 cm⁻¹. These are below the errors with Moss's constant masses (up to 0.28 cm⁻¹).

**Bearing on our claims.** [JK08] tabulates R-dependent mass functions for H₂⁺, but **model (LCAO/FH) functions, not accurate ones**. Its rotational function is off by a factor 1.3–1.9 near R_e, and its level errors are 0.01–0.5 cm⁻¹. Relative to [JK08] our functions are a large quantitative improvement in the same convention. Relative to [S01] they are not. The existence of [S01], which is not cited in [K] or [JK08], explains why [K] still listed "exact wave functions" as future work. It does not make the exact functions new.

## 10. Updated overall novelty statement (supersedes the verdicts in Sections 4 and 7)

Verdict: outcome **(a)**, nearly. Accurate second-order (Bunker–Moss) vibrational and rotational mass functions for H₂⁺ were computed ab initio by Schwenke [S01] in 2001 and tested against Moss's exact non-adiabatic level energies to about 10⁻³ cm⁻¹. That work had a more complete effective Hamiltonian than our first-order mass-only test. Jaquet and Kutzelnigg [JK08] later tabulated less accurate model versions in Kutzelnigg's equivalent convention. Our functions agree with [S01]'s Fig. 1 to the digitization limit (≤ 0.004 everywhere from 1.5 to 19.5 bohr).

* **Not new (now):** the formalism ([BM] 1977, Herman–Asgharian 1966); the qualitative shapes ([K] 2007, [JK08]); accurate ab initio functions for H₂⁺ and HD⁺ over R ≈ 1–20 bohr ([S01], figures plus EPAPS tables); the parameter-free reproduction of the exact H₂⁺ vibration-rotation manifold to 10⁻³ cm⁻¹ ([S01], Table I); the identification of Moss's constants with the functions near R ≈ 3 bohr ([S01]).
* **Ours, as a reproduction with improvements:** an independent recomputation from the exact BO wave function by a different method (Sternheimer equation, spheroidal basis, 3×10⁻¹¹ convergence, 0.02–60 bohr, analytic R⁻³ tail and small-R power laws). [S01] used a Gaussian basis, finite-difference derivatives and a discrete sum over states, with no stated error bar. This is an independent check that confirms [S01] at the level visible in his figure, and it provides a precise benchmark. Minor additions: the first-order mass-only test over all 421 entries of Moss's Table 4 (N ≤ 40), which shows that the mass terms alone carry the non-adiabatic corrections at 10⁻³ cm⁻¹ in our partition; two probable misprints in [M96] Table 4; the reproduction of [JK08]'s deviations, which locates their errors in the model functions; the agreement with m_v(R_e), m_r(R_e) inferred from exact three-body energies at other masses.
* **New (as far as these six papers go):** (1) the observation that in the μ_e-scaled convention one pair of functions serves the entire {X⁺,X⁺,e⁻} family, and its use, with no fitted parameter, to explain the near-constant BO offsets of the threshold critical masses (0.528, 0.512, 0.502 against 0.5375, 0.5062, 0.4970) and the κ⁶ coefficient; (2) the link between the hump of m_v(R) at 3–8 bohr and the effective m_v ≈ 0.49 required by the threshold offsets. Nothing in [BM], [BMM], [K], [M96], [S01] or [JK08] concerns critical masses or mass scaling across the family.
* **Still not checked:** the EPAPS file of [S01] (needed for a sharper than graphical comparison, and to see his R grid); later ab initio work after 2008 (e.g. Pachucki–Komasa on H₂, Mátyus and Teufel 2019 on mass-correction functions; from memory, unverified), which may also contain accurate H₂⁺ functions; the cancellation between [S01]'s C_v and C terms in our partition.

## 11. Suggested wording

**Section 6.2, paragraph "Non-adiabatic mass functions".** Replace the sentence "Kutzelnigg [109] estimated them for H$_2^+$ with LCAO wave functions and wrote down, but did not solve, the exact defining equation." with:

> Schwenke [113] computed them ab initio for H$_2^+$ and HD$^+$ with Gaussian basis sets and showed that the resulting effective Hamiltonian reproduces the exact non-adiabatic H$_2^+$ levels of Moss [112] to about $10^{-3}$ cm$^{-1}$; Kutzelnigg [109] and Jaquet and Kutzelnigg [114] derived an equivalent single-surface formulation and tabulated model (LCAO) versions of the functions.

Replace "We evaluated these functions from the exact BO wave function by solving that equation (...)" with:

> We recomputed these functions independently from the exact BO wave function by solving the defining Sternheimer equation in a prolate-spheroidal basis (2597 values of $R$ from 0.02 to 60 bohr, basis convergence $3\times10^{-11}$). They agree with Schwenke's published curves to the precision of his figure (within 0.004 from 1.5 to 19.5 bohr).

Replace "These shapes agree with Kutzelnigg's LCAO estimates, ... (0.31 and 0.14)." with:

> The model functions of Jaquet and Kutzelnigg [114] have the same shapes but a vibrational maximum of 0.60 at 5 bohr and, at $R_e$, rotational values 1.3 to 1.9 times ours; inserted in the same first-order formula they reproduce their own level deviations (rms 0.003 cm$^{-1}$), so their errors of 0.01 to 0.5 cm$^{-1}$ come from the model functions.

Replace from "Moss [112] fitted constant effective ..." to the end of the paragraph with:

> Moss [112] fitted constant effective vibrational and rotational masses to the exact non-adiabatic level energies of H$_2^+$ and could not determine an $R$ dependence from the energies; Schwenke [113] showed that the ab initio functions supply it. As a check of our functions, inserted to first order with BO vibrational functions they reproduce 421 of the 423 bound entries of Moss's table of non-adiabatic corrections ($v\le19$, $N\le40$) to within 0.004 cm$^{-1}$ (rms 0.001 cm$^{-1}$), comparable to Schwenke's full effective-Hamiltonian result (rms 0.0006 cm$^{-1}$ on dissociation energies); the two exceptions are probably misprints. The functions themselves are therefore not new; what is new here is the observation that, in the $\mu_e$-scaled convention, one pair of functions serves the whole family, and their use, with no fitted parameter, to account for the threshold offsets and the $\kappa^6$ coefficient.

**Abstract clause.** The current wording ("...by the second-order (Bunker–Moss) non-adiabatic vibrational and rotational mass functions, which we evaluate from the exact Born–Oppenheimer wave function.") is acceptable only if it does not imply priority. Safer:

> ...by the known second-order (Bunker-Moss) non-adiabatic vibrational and rotational mass functions of H$_2^+$, which we recompute from the exact Born-Oppenheimer wave function and which, in mass-scaled form, apply to the whole family.

**Introduction line 77** ("a comparison of our mass functions with published ones" outstanding) can be dropped or changed to "a comparison with Schwenke's tabulated functions (EPAPS) beyond graphical accuracy".

**References as printed in the PDFs** (numbers 113 and 114 are suggestions):

* 113. D. W. Schwenke, "A first principle effective Hamiltonian for including nonadiabatic effects for H$_2^+$ and HD$^+$", J. Chem. Phys. **114**, 1693–1699 (2001). doi:10.1063/1.1334897. (Supplementary tables: EPAPS Document No. E-JCPSA6-114-005105.)
* 114. R. Jaquet, W. Kutzelnigg, "Non-adiabatic theory in terms of a single potential energy surface. The vibration–rotation levels of H$_2^+$ and D$_2^+$", Chem. Phys. **346**, 69–76 (2008). doi:10.1016/j.chemphys.2008.02.068.

The note on ref. 109 ("comparison with our $m_v(R)$, $m_r(R)$ still to be made") is out of date: the comparison is in Sections 3, 8 and 9 here.


## 12. Schwenke 2001 EPAPS tables (2026-11-15)

The supplementary file (EPAPS E-JCPSA6-114-005105, available from the journal; not redistributed here, place it at `literature_data/schwenke2001_EPAPS_h2p_tab.txt` to rerun the comparison) tabulates C_vv, C_v, C, C_rr for H2+ (aug-cc-pV5Z) at 33 R values from 0.5 to 200 bohr, and the same for HD+. `schwenke_epaps_compare.py` converts with m_v = -C_vv M^2, m_r = -C_rr M^2 R^2 (M = 1836.15267; the M he used is not printed) and compares with results_fine.csv (output `schwenke_epaps_compare.txt`):

* agreement over 0.5-20 bohr: max |diff| 0.0025 (m_v, near the maximum at 6-7 bohr) and 0.0014 (m_r); median relative difference 0.16% (m_v) and 0.23% (m_r). Schwenke's m_v maximum 0.6581 at 6.0 bohr (ours 0.6603).
* asymptote: Schwenke's functions level off at 0.4989 (both m_v and m_r, R = 100-200), 0.2% below the exact 1/2 (reduced mass of X+ and the atom); ours reach 0.50008 at R = 100 (decaying to 1/2 as R^-3). The constant 0.2% deficit at large R, where the exact value is known analytically, is consistent with a basis-set limit of his Gaussian sum-over-states (or with a slightly different M; a 0.1% smaller M would be needed, which is implausible). Our functions are converged to 3e-11 and reproduce the exact limit.
* Verdict unchanged (Section 10): Schwenke has priority for accurate functions; ours are an independent recomputation that agrees to about 0.2% and is more precise (exact asymptote).

## 13. Mátyus 2018 [EM18] (2026-11-15)

[EM18] E. Mátyus, "Non-adiabatic mass correction to the rovibrational states of molecules: Numerical application for the H$_2^+$ molecular ion", J. Chem. Phys. **149**, 194111 (2018), doi:10.1063/1.5050401 (arXiv:1811.01186). Read in full. Companion: Paper II, J. Chem. Phys. **149**, 194112 (2018), on He$_2^+$ (not read). Numbers below from `matyus2018_compare.py`.

**What it computes.** Most of the paper (Secs. II, III, V) is formal: the second-order non-adiabatic Hamiltonian of Panati–Spohn–Teufel (laboratory-frame Cartesian mass-correction tensor A_{ia,jb} = 2⟨∂ψ|(H−E)^{-1}(1−P)|∂ψ⟩, Eq. 5) is transformed to general curvilinear coordinates for use in the GENIUSH rovibrational program. For a homonuclear diatomic it reduces to δm_vib = A^ρ_ρ = B − b and δm_rot = A^Ω_Ω = A − a, entering as m_n → m_n + δm per nucleus (Eqs. 83–86). Numerics (Sec. IV, VI): floating ECGs with stochastic variational optimization (in-house QUANTEN code, electronic energy converged to < 0.1 μE_h); the resolvent is applied by a linear solve in an auxiliary ECG basis (Σg and Πu) with a complex shift ε = i·10⁻⁶ (a Sternheimer-type solve, not a sum over states); wave-function derivatives by finite differences (Cencek–Kutzelnigg rescaling). R range 0.1–100 bohr. Stated accuracy: non-adiabatic level corrections converged to < 0.005 cm⁻¹; a 10–20 μE_h error in ψ changes the mass functions by < 0.5%; the large-R vibrational function needs a large auxiliary basis. No error bar on the functions themselves.

**What she shows.** (i) The full Cartesian tensor at ρ = 3 bohr (Eq. 67) and the resulting δm_vib = 0.424, δm_rot = 0.144 (Eqs. 84, 86). (ii) Fig. 1: both curves from 0 to 100 bohr with the line 1/2 m_e ("equal share of the single electron between the two atomic nuclei"); vibrational maximum about 0.66 near 5–6 bohr. (iii) Values deposited in her supplementary material (not available here; no table in the paper). (iv) Table I: deviations from Moss's exact non-adiabatic corrections [M96, Moss 1993] for all bound and long-lived quasibound levels, v ≤ 19, J ≤ 41, from solving the full radial equation (Eq. 90, DVR, with DBOC and the ∂A/∂ρ term; m_p/m_e = 1836.152701). Nearly all deviations are −0.001 to −0.004 cm⁻¹. She flags (v, J) = (0, 2) and (10, 4) as unexplained, "might be due to a typographical error".

**Convention.** δm̃_vib = m_v and δm̃_rot = m_r exactly (per nucleus, M → M + m, electron mass 1). The Cartesian elements obey A + a = B + b = 1/2 at ρ = 3 (translation carries the whole electron), a consistency check we do not have in our partition.

**Numerical comparison** (ρ = R = 3 bohr; she prints three decimals):

| | [EM18] | ours | diff | [S01] EPAPS |
|---|---|---|---|---|
| m_v | 0.424 | 0.42538 | −0.0014 (−0.33%) | 0.42552 |
| m_r | 0.144 | 0.14356 | +0.0004 (within rounding) | 0.14346 |

Her m_r agrees with ours to her printed precision. Her m_v is 0.3% low at 3 bohr, outside rounding, while Schwenke and we agree to 0.03% there; this is consistent with her remark that the vibrational function is the harder one to converge in the auxiliary basis. Fig. 1 matches ours by eye (our maximum 0.6606 at 5.91 bohr; m_r = 0.4865, 0.4973, 0.4989 at 10, 15, 20 bohr); a 0.2% asymptotic deficit like Schwenke's cannot be seen at the resolution of the figure. Level tests: her deviations for (0,2) and (10,4), +0.059 and −0.273 cm⁻¹ (Moss minus computed), match ours (+0.061, −0.271), so the two "misprints" we found in Section 7 were identified first by her (she calls them unexplained, possibly typographical). For the rest, her deviations are systematically −0.002 to −0.003 cm⁻¹ (her values larger than Moss's), ours are +0.0007 mean, rms 0.001 (ours slightly smaller); the 0.002–0.003 cm⁻¹ difference between the two treatments is of the size of the extra O(M⁻²) pieces she includes and we omit (full equation with ∂A/∂ρ term and DBOC-modified χ, versus our first-order expectation value), and within her stated 0.005 cm⁻¹ accuracy.

**Thresholds, families, large-R limit.** She notes the large-R limit 1/2 as the equal-share value (one sentence in the Fig. 1 caption) and, from the translational block, that the total mass becomes the molecular mass (Eq. 81). She does not discuss mass scaling, the {X⁺,X⁺,e⁻} family, the μ_e-scaled argument, critical masses or thresholds. She proposes comparing with Korobov's three-body resonances as future work (Sec. VII).

**Novelty verdict (update of Section 10).** [EM18] pre-empts more than [S01] in three respects: (1) an ab initio Sternheimer-type (resolvent by linear solve, not sum over states) evaluation of the H$_2^+$ functions over 0.1–100 bohr, i.e. the method class we use; (2) a full level-by-level test against Moss's table for all v ≤ 19, J ≤ 41, including quasibound levels, at a few 10⁻³ cm⁻¹, which is the same test as ours (ours is first order and covers only bound levels); (3) the identification of the two anomalous Moss entries. Our "minor additions" in Section 10 must therefore drop "the first-order mass-only test over all of Moss's Table 4" and "two probable misprints" as new items; at most we confirm them. What remains ours: an independent evaluation in a spheroidal basis with stated convergence (3×10⁻¹¹) and the exact asymptote, more precise than both [S01] (0.2% asymptotic deficit) and [EM18] (m_v 0.3% low at 3 bohr, functions not tabulated in the paper); and, unchanged, items (1) and (2) of Section 10 "New" (one pair of functions for the whole family in μ_e-scaled form; parameter-free account of the threshold offsets and the κ⁶ coefficient; the hump of m_v linked to the effective 0.49). Neither [S01] nor [EM18] touches critical masses or mass scaling. Still not checked: her supplementary values (would allow a 0.1% comparison over R), Paper II, Pachucki–Komasa 2009 (H2, same functions), Scherrer et al. PRX 2017.

**Suggested wording.**

Section 6.2, replacing the sentence that credits Schwenke [113] (Section 11 above):

> Schwenke [113] computed them ab initio for H$_2^+$ and HD$^+$ with Gaussian basis sets, and Mátyus [115] recomputed them for H$_2^+$ with floating explicitly correlated Gaussians from 0.1 to 100 bohr; both showed that the resulting effective Hamiltonian reproduces the exact non-adiabatic H$_2^+$ levels of Moss [112] to a few $10^{-3}$ cm$^{-1}$, and Mátyus noted that two entries of Moss's table are probably misprinted. Kutzelnigg [109] and Jaquet and Kutzelnigg [114] derived an equivalent single-surface formulation and tabulated model (LCAO) versions of the functions.

And, for the agreement sentence: "They agree with Schwenke's tabulated values to about 0.2% and with the values printed by Mátyus at 3 bohr to 0.3% or better, and they reach the exact separated-atom limit of one half."

In the Moss paragraph (Section 11), replace "the two exceptions are probably misprints" by "the two exceptions are the entries that Mátyus [115] already identified as probably misprinted", and "comparable to Schwenke's full effective-Hamiltonian result" by "comparable to the full effective-Hamiltonian results of Schwenke [113] and Mátyus [115]".

Abstract clause:

> ...by the known second-order (Bunker-Moss) non-adiabatic vibrational and rotational mass functions of H$_2^+$ [113,115], which we recompute from the exact Born-Oppenheimer wave function and which, in mass-scaled form, apply to the whole family.

Reference as printed (number 115 is a suggestion):

* 115. E. Mátyus, "Non-adiabatic mass correction to the rovibrational states of molecules: Numerical application for the H$_2^+$ molecular ion", J. Chem. Phys. **149**, 194111 (2018). doi:10.1063/1.5050401. (Supplementary material: mass-correction values for H$_2^+$.)
