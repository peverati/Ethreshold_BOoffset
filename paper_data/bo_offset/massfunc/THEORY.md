# Nonadiabatic mass functions m_v(R) and m_r(R) for {X+, X+, e-}: theory
*Last updated: 2026-11-15 19:02*

Written 2026-11-14. Derivation of the formulas coded in `massfunc.py`. Everything not marked otherwise is derived here; literature statements are flagged "(memory, unverified)".

## 1. Hamiltonian and conventions (same as solvers/bo_h2plus)

Jacobi coordinates: **R** = X–X vector, **ρ** = electron position relative to the X–X midpoint (the nuclear centre of mass; no mass-polarization term). Nuclear mass M, electron mass 1:

    H = T_N + H_el(μ_e),   T_N = −(1/2μ_N) ∇_R²,   μ_N = M/2,
    H_el(μ_e) = −(1/2μ_e) ∇_ρ² + 1/R − 1/|ρ − R/2| − 1/|ρ + R/2|,   μ_e = 2M/(2M+1).

BO in the project convention uses H_el(μ_e). With φ0 the 1sσg eigenfunction (energy E0(R)) and the unit-mass problem h, e0(r), φ(ρ; r):

    φ0(ρ; R) = μ_e^{3/2} φ(μ_e ρ; μ_e R),   E0(R) = μ_e e0(μ_e R).            (1)

The zero-energy (threshold) radial problem of `paper_data/bo_offset/README.md` is written in r = μ_e R with coupling λ = M/μ_e = M + 1/2. "M → M + m" there means λ → λ + m in the kinetic term.

## 2. Second-order effective nuclear Hamiltonian

Partition Ψ = φ0 χ + (Q part), P = |φ0⟩⟨φ0| acting on ρ at each **R**. Löwdin partitioning gives

    H_eff = P H P + P T_N Q [E − Q H Q]^{-1} Q T_N P.

P H P is BO + diagonal (adiabatic) correction. For the coupling,

    Q T_N (φ0 χ) = −(1/μ_N) Q(∂_i φ0) ∂_i χ + Q(T_N φ0) χ.

Take the leading resolvent, [E − QHQ]^{-1} ≈ −Q (H_el − E0)^{-1} Q. Its corrections are relatively O(1/M) because nuclear kinetic and vibrational energies are O(M^{-1/2}) and smaller than electronic gaps. Keep the term with two nuclear derivatives:

    ΔH = −(1/μ_N²) ∂_i† A_ij ∂_j = (1/μ_N²) ∂_i A_ij(R) ∂_j,
    A_ij(R) = ⟨∂_i φ0 | Q (H_el − E0)^{-1} Q | ∂_j φ0⟩_ρ  (∂ at fixed ρ).                 (2)

Combine it with T_N = −(1/2) ∂_i (1/μ_N) ∂_i:

    T_N + ΔH = −(1/2) ∂_i [ δ_ij/μ_N − 2A_ij/μ_N² ] ∂_j = −(1/2) ∂_i [(μ_N + 2A)^{-1}]_ij ∂_j + O(M^{-3}).

So **the nuclear reduced mass becomes the tensor μ_N 1 + 2A(R)**. A is positive semidefinite because H_el − E0 > 0 on Q. So the mass always increases, binding increases, and the exact critical masses lie *below* BO+ad. This matches the sign of the observed offsets.

Neglected (all O(M^{-2}) in the energy, i.e. O(1/M) in λ_c, so they fall off with v):

* the cross and diagonal second-order terms with Q(T_N φ0) (potential-like);
* the O(1/M) corrections to the resolvent;
* third-order terms.

## 3. Symmetry: parallel (Σg) and perpendicular (Πg) parts

Cylindrical symmetry gives A_ij = A_∥ n_i n_j + A_⊥ (δ_ij − n_i n_j), with n = R/R. Take R along z.

* **Parallel:** f_∥ = ∂_R φ0|_ρ. It is σ (m = 0). It is **gerade**, because inversion ρ → −ρ commutes with ∂_R at fixed ρ. Only Σg states couple. The 2pσu state, which becomes degenerate with 1sσg as R → ∞, is excluded by symmetry. Without that, the resolvent would diverge.
* **Perpendicular:** ∂_x φ0|_ρ comes from rotating the frame. It equals −(1/R)(z∂_x − x∂_z)φ0 = ±(i/R) L_y φ0 (electron angular momentum about the midpoint). For axial φ0(s, z): (z∂_x − x∂_z)φ0 = −cos φ · g, with g = s ∂_z φ0 − z ∂_s φ0. In spheroidal coordinates

      g = sqrt((ξ²−1)(1−η²)) (η φ_ξ − ξ φ_η)/(ξ² − η²).

  L_y is even under inversion, so f_⊥ = g cos φ / R is **Πg**. Only Πg states couple (no projector needed).
* **Check:** ⟨f_∥|f_∥⟩ + 2⟨f_⊥|f_⊥⟩ = D(R), the bo_h2plus adiabatic function. `massfunc.reference_norms` and the `fpar_norm_in_basis` / `fperp_norm_in_basis` columns test this.

## 4. Mass scaling and the radial equation in the project convention

Use (1): ∂_R → μ_e ∂_r and H_el − E0 → μ_e (h − e0). Then A^{(μe)}(R) = μ_e A(μ_e R), where A is the unit-mass quantity. With nuclear angular momentum L, the radial equation of T_N + ΔH + V is

    −(1/2) d/dR [1/(μ_N + 2A_∥)] d/dR u + L(L+1)/(2(μ_N + 2A_⊥)R²) u + V u = E u.

Multiply by M/μ_e² and set r = μ_e R. The radial kinetic prefactor becomes

    M/(M + 4 μ_e A_∥(r)) = λ/(λ + 4 A_∥(r))   (exact algebra, λ = M/μ_e),

and the same holds for the centrifugal term with A_⊥. Hence, with no further approximation inside the effective Hamiltonian,

    **m_v(r) = 4 A_∥(r),   m_r(r) = 4 A_⊥(r)**   (unit-mass A, r = μ_e R),   λ → λ + m.              (3)

Equivalently, in kappa6.md's convention (−(1/M) d²/dR², potential μ_e U(μ_e R)), M → M + 4μ_e A_∥(μ_e R) ≈ M + m_v(R) + O(1/M).

## 5. Separated-atom limit m_v(∞) = m_r(∞) = 1/2 (derived)

For R → ∞, φ0 → (a + b)/√2, with b = 1s centred at +R/2 ẑ.

* **Parallel:** ∂_R b|_ρ = −(1/2) ∂_z b, so f_∥ → −(1/2)(∂_z a' + ∂_z b)/√2, with the two atomic pieces non-overlapping.
* **Perpendicular:** L about the midpoint is L_B + (R/2) ẑ × p. Since L_B 1s = 0, f_⊥ → −(1/2) ∂_y-type atomic pieces as well.
* **Both reduce to the atomic quantity** (1/4)⟨∂_z 1s|(h − e)^{-1}|∂_z 1s⟩. With p_z = i[h, z] (unit mass), (h − e)^{-1} p_z |1s⟩ = i z|1s⟩ (it is automatically orthogonal to 1s). So

      ⟨p_z 1s|(h − e)^{-1}|p_z 1s⟩ = i⟨1s|p_z z|1s⟩ = i(−i/2) = 1/2,

  and A_∥, A_⊥ → 1/8. Then **m_v(∞) = m_r(∞) = 4/8 = 1/2**.
* **Physical check (exact):** the X+ + (Xe) relative motion has 2μ = 2M(M+1)/(2M+1), so λ_∞ = 2μ/μ_e = M + 1 = (M + 1/2) + 1/2. The asymptotic shift 1/2 is exact in this convention.
* Corrected 2026-11-15: in a clamped-nuclei convention (Kutzelnigg 2007, 1/μ_vib = (2/M)(1 + A/M), m_v = −A) the asymptote is also m(∞) = 1/2, since the exact asymptotic reduced mass M(M+1)/(2M+1) corresponds to M + 1/2 per nucleus. The earlier statement that it doubles to 1 was wrong (a value of 1 arises only for two electrons, H2, or per total nuclear mass). See LITERATURE_COMPARISON.md.

## 6. Numerical method (Sternheimer / Hylleraas, no sum over physical states)

A = ⟨f|Q (h − e0)^{-1} Q|f⟩ = −min_ψ J[ψ], with J[ψ] = ⟨ψ|h − e0|ψ⟩ + 2⟨ψ|Q f⟩ (Hylleraas). The minimizer solves the Sternheimer equation (h − e0)ψ1 = −Q f. In a finite basis the result is a **lower bound** that converges from below.

* **Basis** (prolate spheroidal ξ, η, φ; x = 2α(ξ − 1), α = p = R sqrt(−e0/2)):
  * Σg: L_i(x) e^{−x/2} P_l(η), l even.
  * Πg: sqrt(ξ²−1) L_i(x) e^{−x/2} P_l^1(η) cos φ, l even (gerade needs h odd in η for m = 1).
* **Integrals are exact.** The kinetic, overlap, Coulomb and m² centrifugal integrands are polynomial × e^{−x}, using dV = (R/2)³(ξ²−η²), 1/r_A + 1/r_B = 4ξ/(R(ξ²−η²)) and (ξ²−η²)/((ξ²−1)(1−η²)) = 1/(1−η²) + 1/(ξ²−1). So Gauss–Laguerre × Gauss–Legendre quadrature is exact. Every matrix is a sum of Kronecker products of 1D matrices.
* **Solve.** The generalized eigenproblem H c = ε S c is solved in the basis. Then A = Σ_{n≥1} ⟨n|f⟩²/(ε_n − ε_0), with n = basis pseudostates. This is identical to solving the Sternheimer linear system in the basis, so it is not a slowly converging sum over physical states. The basis ground state (ε_0 − e0 ≤ 1e−12 for R ≥ 0.5) is removed (Q).
* **Right-hand sides.** f_∥ and g come from the exact bo_h2plus wavefunction (Jaffé × Legendre, about 1e−12), with the chain rule of `bo_h2plus/adiabatic.py`:

      f_∥ = ∂_R φ|_{ξη} − [ξ(ξ²−1)φ_ξ + η(1−η²)φ_η]/(R(ξ²−η²)).

  ∂_R|_{ξη} is a 5-point difference, h = 1e−3 R. The (ξ²−η²) volume factor cancels the 1/r_{A,B} kinks, so the projections ⟨χ_k|f⟩ have smooth integrands.

## 7. Using m(r) for the critical-mass offsets

**First order** (as in the parent README). With u the BO+ad zero-energy solution at λ_c and B = ∫(−ΔU)u²,

    M_c^{BO+ad} − M_c^{exact} ≈ [∫ m_v(r) u'² dr + ∫ m_r(r) L(L+1) u²/r² dr] / (λ_c B).                 (4)

**Nonperturbative variant** (sensitivity check, `analyze.py`). Solve

    −(q_v u')' + q_r L(L+1)/r² u + λ ΔU u + ΔD u = 0,   q = λ/(λ + m(r)),

for λ_c directly. Its difference from (4) is O(m²/λ). That is formally the same order as the neglected terms of Section 2, so the spread between the two is an honest error bar for the "no-fit" prediction.

## 8. Uncertain points

* **Sign of g and of the perpendicular derivative:** irrelevant, because A is quadratic in f.
* **Normalization of P_l^1 and the cos φ factor:** the φ integral (π for cos²φ, 2π for m = 0) multiplies b^T M^{-1} b once. Checked numerically through ⟨f|f⟩ against D(R) (agreement to 1e−4 relative, limited by basis completeness for f, not for A).
* **Identification with literature "vibrational/rotational masses"** (Bunker–Moss; Kutzelnigg Mol. Phys. 2007; Pachucki–Komasa JCP 2009): conventions (electron from nuclear CM versus geometric centre, μ_e versus 1, and whether the 1/2 is counted) differ. Not checked against their tables (memory, unverified). For H2+ these coincide (midpoint = nuclear CM).
* **Neglected terms** (Section 2) are O(1/M) in λ_c. They matter most for L = 1 (M_c ≈ 4.2) and v = 1 (M_c ≈ 9.5).
