# Critical masses of {A+, B+, e-} with unequal nuclear masses: theory

*Last updated: 2026-11-16 10:11*

Written 2026-11-16. Everything here is derived in this file and checked numerically where stated (tests in
`solvers/bo_h2plus/tests/test_unequal.py` and `test_unequal_masses.py`). Literature statements are flagged.

## 1. System, coordinates, mass constants

Nuclei A (mass m_A = M, the light one) and B (m_B = r M, r ≥ 1), electron mass 1, charges +1, +1, −1.
Jacobi coordinates: **R** = R_B − R_A, **ρ** = electron relative to the nuclear centre of mass. Then

    H = −(1/2μ_N) ∇_R² − (1/2μ_e) ∇_ρ² + 1/R − 1/r_A − 1/r_B,
    μ_N = rM/(1+r),   M_eff ≡ 2μ_N = 2rM/(1+r),   μ_e = (1+r)M/((1+r)M + 1).

There is no mass-polarization term. The nuclear centre of mass lies at **+βR** from the midpoint (towards B),

    β = (m_B − m_A)/(2(m_A + m_B)) = (r − 1)/(2(1 + r))       (r = 2: β = 1/6; r = 3/2: β = 1/10).

**Thresholds (exact).** The two atoms have reduced masses μ_Ae = M/(M+1) < μ_Be = rM/(rM+1). The lowest
continuum (for every natural-parity L) is A⁺ + (B e)(1s), **E_th = −μ_Be/2**; A⁺... the other channel,
B⁺ + (A e)(1s), lies (μ_Be − μ_Ae)/2 ≈ (r−1)/(2rM) higher. The open-channel ion–atom reduced mass is
μ_rel = M(rM+1)/((1+r)M+1).

## 2. The zero-energy (coupling-constant) form survives, with a new linear map

Clamped-nuclei problem: the electronic Hamiltonian at fixed R is that of H2+ (the nuclear charges are equal),
in midpoint coordinates ρ_m = ρ + βR. With electron mass μ_e, φ(ρ_m; R) = μ_e^{3/2} φ₁(μ_e ρ_m; μ_e R) and
E(R) = μ_e U(μ_e R) exactly (U = unit-mass H2+ curve incl. 1/R). Write the radial equation in r = μ_e R and
multiply by M_eff/μ_e²:

    −u'' + L(L+1)/r² u + λ U(r) u + [O(1) terms] u = (M_eff/μ_e²) E u,      λ ≡ M_eff/μ_e.

**λ is again linear in M:** λ = 2rM/(1+r) + 2r/(1+r)², i.e.

    **M_c = (1+r) Λ_c/(2r) − 1/(1+r)**          (r = 1: M_c = Λ_c − 1/2, the symmetric result).

So the critical mass is still an eigenvalue of an M-independent linear problem; any O(1) term in the bracket
(order 1/M in energy) shifts Λ_c by an M-independent amount, and the shift in M is multiplied by (1+r)/(2r)
(= 3/4 for r = 2). The **BO critical couplings Λ_c^BO are identical to those of the symmetric family**
(same U), only the map to M changes: M_c^BO(r) = (1+r)Λ_c^BO/(2r) − 1/(1+r).

## 3. Diagonal (adiabatic) correction of 1sσg

At fixed ρ the nuclear derivatives acquire an electronic derivative (R along z, B at z = +R/2):

    ∂_R|_ρ = ∂_R|_{ρ_m} + β ∂_z,        ∂_{R_x}|_ρ = (1/R)(x∂_z − z∂_x) + β ∂_x.

Cross terms between the gerade (first) and ungerade (second) pieces vanish in diagonal elements, so

    K_gg = ⟨∇_R φ_g|∇_R φ_g⟩_ρ = D(r) + β² t_g(r),     t_g = ⟨p²⟩_g = −2(U + rU')  (virial),

with D the bo_h2plus function (D(∞) = 1/4) and t_g(∞) = 1. **BO+ad (g only, own asymptote)**:

    −u'' + L(L+1)/r² u + λ ΔU u + [ΔD + β²(t_g − 1)] u = 0,    ΔU = U + 1/2, ΔD = D − 1/4.

Check (exact algebra): in clamped-nuclei units (electron mass 1, electron from the geometric centre) the same
O(1/M) diagonal correction is (1/m_A + 1/m_B)(D/2 + t_g/8): the usual DBOC of HD⁺ is the H2⁺ one times
(1/m_p + 1/m_d)/(2/m_p). (Compared with Schwenke 2001's HD⁺ adiabatic table, `literature_data`, in the tests.)

**The asymptote of this g-only model is wrong for r ≠ 1.** Its limit −λ/2 + 1/4 + β² equals the *mean* of the
two atomic thresholds to O(1/M). In scaled units the exact open threshold lies **β below** it and the closed one
β above (exact values: `Masses.thresholds_scaled`). Two g-only conventions are therefore reported:

* **BOad** — measured from its own asymptote (continuous with the symmetric family at r = 1);
* **BOad_thr** — measured from the exact threshold: the g level must be bound by β (scaled) below the g
  asymptote. For L = 0 this removes the zero-energy halo (the g channel is closed, κ² = β).

## 4. The g/u near-degeneracy: two-channel model (non-perturbative)

The coupling β∂_z connects 1sσg to σu states. 2pσu becomes degenerate with 1sσg as R → ∞, so second-order
perturbation theory fails exactly where the asymmetric threshold is decided. 2pσu is therefore treated
**explicitly**: Ψ = φ_g χ_g + φ_u χ_u. With a_ij = ⟨φ_i|∂_R φ_j⟩_ρ (radial; antisymmetric) and
K_ij = ⟨∇_R φ_i|∇_R φ_j⟩_ρ, the projected nuclear kinetic energy gives (R̂-radial derivative coupling only;
the perpendicular ⟨φ_g|∂_{R_x} φ_u⟩ vanishes by cylindrical symmetry):

    −u'' − 2a u' − a' u + [L(L+1)/r² + K(r)] u + λ diag(U_g, U_u) u = (M_eff/μ_e²) E u,

with, by parity,

    a_gu = β ⟨g|∂_z|u⟩ = β (E_u − E_g)⟨g|z|u⟩ (hypervirial; checked to 1e-9),
    K_gg = D_g + β² t_g,  K_uu = D_u + β² t_u,
    K_gu = β k_gu,  k_gu = ⟨F_g|Z_u⟩ + ⟨Z_g|F_u⟩ + 2(⟨G_g|S_u⟩ + ⟨S_g|G_u⟩)/R,

F = ∂_Rφ|_{ρ_m}, Z = ∂_zφ, S = ∂_sφ, G = s∂_zφ − z∂_sφ (all from the exact bo_h2plus wave functions;
`bo_h2plus/unequal.py`). Limits (derived, checked numerically): D_u → 1/4, t → 1, k_gu → −1, a_gu → 0.
Hence K(∞) = [[1/4+β², −β], [−β, 1/4+β²]], eigenvalues (1/2 − β)² = 1/(1+r)² (eigenvector (g+u)/√2 =
electron on B) and (1/2 + β)² = r²/(1+r)² (on A): exactly the two atomic reduced-mass shifts. **The two-channel
model has the exact threshold built in** (open channel at −β, closed at +β relative to the g asymptote, to O(1)),
and at r = 1 it decouples into the symmetric BO+ad problem.

All terms are λ-independent or ∝ λ, so the critical couplings are still eigenvalues of an M-independent
problem and §2 holds. But the g/u mixing is **not** a small perturbation: it is set by the competition between
the exchange splitting λ(U_u − U_g) (exponentially small in r) and β. The electron localizes on B beyond
r_loc(λ) where λ(U_u − U_g) ≈ 2β; r_loc grows like ln λ. Consequences (predictions):

1. Λ_c(2ch) lies between Λ_c(BOad) and Λ_c(BOad_thr), far from both; the offset relative to BOad is large
   (several units of λ) and **grows with v** (the higher the level, the smaller dE/dλ at threshold, the larger
   the λ needed to supply the extra binding β, partly returned by the localization).
2. For β → 0 the shift goes as √β (L = 0: a closed component must decay at the threshold; checked numerically:
   r = 1.001 and 1.01 give shifts in ratio 3.05 ≈ √10).

## 5. Non-adiabatic masses in the {g, u} model space (second order, parameter-free)

Löwdin partitioning with model space {g, u} (Q = 1 − |g⟩⟨g| − |u⟩⟨u|); the leading term with two nuclear
derivatives gives an R-dependent **mass matrix** (as in massfunc/THEORY.md §2, now 2×2):

    A_ij(e) = ⟨f_i|Q(h − e)^{-1}Q|f_j⟩,  A_ij = [A_ij(e_i) + A_ij(e_j)]/2,  m_ij = 4A_ij,
    radial  f_i = F_i + βZ_i (Σ blocks),  perpendicular f_i = G_i/R + βS_i (Π blocks),
    −(d/dr) Q_∥ (d/dr) and L(L+1) Q_⊥/r²,  Q = (1 + m/λ)^{-1} (non-perturbative; first order is 1 − m/λ).

By parity: m_gg = m_v(r) + β² m_gg2(r) (m_v = the symmetric-family function), m_uu analogous,
m_gu = β m_gu1(r). Computed with the Sternheimer code of massfunc (Σg, Σu, Πg, Πu blocks;
`massfunc_gu.py`, table `massfunc_gu.csv`, basis convergence ≤ 2e-12 for r ≥ 0.1).
Exact checks: m⊥_gg2 ≡ 2 (sum rule (h − e_g)^{-1}∂_x g = −x g, Πu complete) and m∥_gg2 + m∥_uu2 ≡ 4
(the 2pσu part moved between rows); atomic limits m_gg = m_uu = 1/2 + 2β², m_gu = −2β, so the open channel has
m = 1/2 − 2β + 2β² = **2/(1+r)²**, exactly the value required by the exact open-channel reduced mass
(2μ_rel/M_eff · λ − λ → 2/(1+r)²). The u-row pieces diverge as r → 0 (united-atom degeneracy of 2pσu, 2pπu,
2sσg); there the u channel is closed by ≥ λ·0.9, and they are frozen below r = 2 (variant: r = 3) and switched
off as (r/2)². Neglected: potential-like second-order terms and resolvent corrections (relative O(1/λ)),
as in the symmetric theory.

## 6. Models used for the predictions (`twochannel.py`)

| kind | content | role |
|---|---|---|
| BO | g, λΔU_g | baseline (Λ_c identical to symmetric family) |
| BOad | + ΔD + β²(t_g − 1), own asymptote | conventional BO+ad |
| BOad_thr | BOad measured from the exact threshold | BO+ad, correct threshold, no mixing |
| BOad+m | BOad + g-row masses | "symmetric theory transferred" |
| 2ch | g, u, a, K, exact threshold | non-perturbative g/u mixing |
| **2ch+m** | 2ch + 2×2 mass matrices | **parameter-free theory** |
| 2ch+m_asym, 2ch+m_rc3 | u-row masses atomic / frozen below r = 3 | error-bar variants |

Numerics: outward regular solutions from r = 1e-3, inward closed-decaying and open-tail solutions from
r_m = max(40, 30/√(2β)), QR renormalization with sign tracking, matching determinant at r = 6. For one channel
the code reproduces paper_data/bo_offset/coupling.py to 1e-10 (BO, BO+ad, L = 0, 1) and massfunc's non-perturbative
offset (0.5281, v = 1).

## 7. Null models and where theory and null model differ

* **Persistence (referee M4):** offset(v) = offset(v − 1) of the same family, offset measured from a fixed
  reference (BOad or BOad_thr). Evaluated by analyze.py after the exact runs.
* **Symmetric transfer:** offset in λ units = the exact symmetric-family offset of the same (L, v).

The theory predicts offsets from BOad of several units in M that grow with v (the symmetric family: 0.50–0.54,
constant to 0.04). Persistence then fails by the predicted level-to-level change (≫ 0.03); the symmetric
transfer fails by the whole g/u term. Measured from BOad_thr the theory offsets are negative (exact below
BOad_thr: localization helps) and also v-dependent. Numbers: PREDICTION_UNEQUAL.md.

## 8. Uncertain points (honest list)

* The two-channel model omits 2pπu/Πg channels explicitly (they enter only through the masses) and all higher
  Σ states; their effect is second order and O(1/λ) relative, as in the symmetric theory, **but** the
  threshold weights of an L = 0 state at the asymmetric threshold extend into the localization region, where
  the mass functions vary; the variants bracket this.
* O(1/λ) neglected terms: in the symmetric family they produced residuals ≤ 0.012 (L = 0) and up to 0.03
  (L = 1, M_c = 4.2). For r = 2 the λ-to-M factor 3/4 reduces them, but the larger g/u structure may enlarge
  them; the error bar adds the symmetric residual of the same (L, v) to the variant spread.
* The exact critical mass refers to the true threshold −μ_Be/2 (B e + A⁺); levels between the two thresholds
  are bound in the three-body sense only if below −μ_Be/2.
