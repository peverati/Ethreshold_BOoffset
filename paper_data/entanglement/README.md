# Electron-nuclear entanglement of the exact {X⁺, X⁺, e⁻} ground state versus M

*Last updated: 2026-11-13 21:06*

Written 2026-11-13. The exact L = 0 ground-state wave functions come from `paper_data/mass_scans/data/L0_M*.npz` (solvers/threebody, perimetric Laguerre basis). The N−4 and fixed-basis re-solves use `paper_data/fidelity/fidcore.py`. The Born–Oppenheimer product states come from `solvers/bo_h2plus`. All three packages were used unchanged; the code in this folder is new. Atomic units, κ = M^(−1/4).

## Regenerate

```bash
cd paper_data/entanglement
export THREEBODY_BLAS_THREADS=3 OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3; PY=python3
$PY validate.py                      # -> validation.txt (~1 min)
$PY run.py exact                     # all 37 scan masses, base quadrature (~10 min)
$PY run.py fine                      # finer quadrature (convergence)
$PY run.py prevN 0.05 0.4 1 3 20 206.768283 1836.15267343 10000   # N-4 basis (convergence)
$PY run.py bo 3 4.5 6 9.5256 12 15 20 29.639 50 100 206.768283 500 1000 1836.15267343 4000 10000
$PY run.py dense 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.4924 0.5 0.525 0.55 0.575 0.6
$PY analyze.py                       # table.csv, fits.txt, fig1..fig4
```
Each mode appends one row per mass to `results_<mode>.csv` and skips masses already there, so you can resume a run by repeating the command. Each mode also saves the spectra to `data/ent_<mode>_M*.npz`.

## Definitions

* **Bipartition.** R = x₁ − x₀ (the X–X vector, "nuclear") versus ρ = x₂ − (x₀+x₁)/2 (the electron relative to the X–X midpoint, "electronic"). Both are mass-independent Jacobi vectors, and the kinetic energy is separable in them. So R|ρ is the natural electron–nuclear split, and it is the split in which BO is a product state (in the body frame). Distances: r₀₁ = R, r₀₂,₁₂² = ρ² + R²/4 ± ρRu, with u = cos∠(R, ρ).
* **Schmidt structure (exact).** For L = 0, Ψ = Σ_l g_l(R, ρ) Σ_m Y_lm(R̂) Y*_lm(ρ̂), where g_l = 2π ∫du Ψ P_l(u). The reduced density matrix of R is therefore block-diagonal in l. Block l has radial kernel ∫ρ²dρ g_l(R,ρ) g_l(R',ρ), and each of its eigenvalues has multiplicity 2l+1. Exchange symmetry (X₀ ↔ X₁: R → −R) leaves only even l.
* **Measures.** S = von Neumann entropy (nats). S_lin = 1 − Tr γ². Weights w_l = probability of partial wave l, which equals the probability that the electron's orbital angular momentum about the midpoint is l (and also that of the X–X vector, since L = 0).
* **Exact split.** S = S_orient + S_rad (an identity, not a model):
  * S_orient = −Σ w_l ln(w_l/(2l+1)): the distribution of l, plus the m-degeneracy.
  * S_rad = Σ_l w_l S_l, where S_l is the entropy of the normalized radial block l.
  * Also reported: S₀, the radial entropy of the s-wave block alone, and p₂, the second Schmidt coefficient of the l = 0 block.
* **Radial-only entanglement.** Take the partial traces over the directions R̂ and ρ̂. These are local operations, so they cannot increase entanglement. The result is the mixed state ρ_rad = Σ_l (2l+1)|g_l⟩⟨g_l| of (|R|, |ρ|). Two quantities are computed for it:
  * its logarithmic negativity E_N, computed in a compressed local basis; it is an entanglement measure of the radial state and a lower bound on what survives the angular trace;
  * the classical mutual information I(|R|;|ρ|) of P(R, ρ), which counts all correlation, classical and quantum, between the two magnitudes.
* **Orientational coherence.** After tracing out the electron, the angular part of γ_R is Σ_l w_l P_l(R̂·R̂′)/(4π). D(Θ) = Σ_l w_l P_l(cos Θ) is the coherence between two orientations of the X–X axis that differ by an angle Θ, normalized so that D(0) = 1. A pure rotor in l = 0 would have D ≡ 1.
* **BO product state.** Ψ_BO = [P_BO(R)/(4πR²)]^½ φ_1sσg(ρ; R): the exact clamped-nuclei electron times the BO v = 0, J = 0 nuclear density. It is evaluated on the same grids.

## Numerics

* **Quadrature.** Composite Gauss–Legendre in R on [⟨R⟩ − 10ΔR, ⟨R⟩ + 14ΔR] (clipped at 0), in ρ on [0, 14⟨r_Xe⟩] (panels uniform in √ρ), and in θ ∈ [0, π/2] (128 nodes; the coalescence kinks sit at the endpoint). lmax = 64.
* **Fine grid.** About 1.5× more points in every direction, and lmax = 96.
* **Spectrum.** Nyström discretization: SVD of A^l_ij = √a_i g_l(R_i, ρ_j) √b_j.
* **Partial-wave tail.** w_l decays as l⁻⁶, set by the electron–X cusps; this is measured, e.g. in the M → ∞ electronic function between l = 80 and l = 120. The tail beyond lmax is estimated in `table.csv` (`tail_S`, ≤ 1e−5).

## Validation (`validation.txt`)

| test | result |
|---|---|
| separable Gaussian e^(−0.6R²−0.9ρ²), separable exponential e^(−R−2ρ) | S = 1e−26, w₀ = 1, E_N = 0, MI = 1e−16 |
| coupled Gaussian e^(−aR²−bρ²−cR·ρ) (a, b, c = 0.5, 0.7, 0.6): analytic S = 3 S_1D (three Cartesian two-mode Gaussians, Mehler kernel) | S = 0.856460302494, error 5e−15; purity error 2e−15. The 1D closed form was itself checked against a direct 1D grid SVD (agreement to 14 digits). |
| Parseval: Σ_l (2l+1)‖A^l‖² against the direct 6D quadrature norm | agree to the l-tail (≤ 4e−7). The direct norm equals 1 to ≤ 8e−7 at every M. |
| quadrature base → fine (37 masses) | S ≤ 1.2e−5, S_orient ≤ 1.0e−5, S_rad ≤ 2e−6, MI ≤ 2e−6; **E_N ≤ 5.7e−3 for M < 1** (≤ 1e−3 above). This error comes from the compression cut, so E_N is good only to about 0.5% at small M and shows grid jitter. |
| basis N → N−4 (8 masses, 0.05 to 10⁴) | S ≤ 2.5e−7, S_rad ≤ 1.3e−7, MI ≤ 4e−8, E_N ≤ 9e−4 |
| the null test | the solvable harmonic models of `solvers/harmonic_model` are exact products in R\|ρ (S = 0); the separable tests above play this role on the actual code path |
| M → ∞ orientational limit | the BO product state and the exact states approach the same S_orient(∞) (next section) |

**Error budget.** S is converged to about 1e−5 absolute everywhere (relative 1e−5 to 4e−5). The l-tail is of the same size. The basis error is negligible.

## Results

| M | S | S_lin | S_orient | S_rad | w₀ | w₂ | D(90°) | E_N (radial) | I(\|R\|;\|ρ\|) | S (BO product) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 1.495240 | 0.469641 | 1.023761 | 0.471480 | 0.7975 | 0.1545 | 0.7303 | 0.941 | 0.374467 | – |
| 0.4 | 1.603937 | 0.500446 | 1.107124 | 0.496813 | 0.7777 | 0.1677 | 0.7052 | 0.975 | 0.403208 | – |
| 1 (Ps⁻) | 1.516233 | 0.481498 | 1.051035 | 0.465198 | 0.7883 | 0.1637 | 0.7169 | 0.928 | 0.363622 | – |
| 3 | 1.173831 | 0.389987 | 0.819000 | 0.354831 | 0.8380 | 0.1343 | 0.7776 | 0.751 | 0.239628 | 2.4213 |
| 9.5256 | 0.808871 | 0.275414 | 0.575451 | 0.233420 | 0.8910 | 0.0963 | 0.8463 | 0.548 | 0.128245 | 1.0853 |
| 20 | 0.645103 | 0.218911 | 0.469950 | 0.175153 | 0.9136 | 0.0783 | 0.8768 | 0.444 | 0.085701 | 0.7577 |
| 206.77 (Mu₂⁺-like) | 0.383418 | 0.124983 | 0.313876 | 0.069543 | 0.9459 | 0.0507 | 0.9216 | 0.236 | 0.025906 | 0.3935 |
| 1836.15 (H₂⁺) | 0.298770 | 0.095670 | 0.270449 | 0.028321 | 0.9545 | 0.0430 | 0.9338 | 0.133 | 0.008711 | 0.3008 |
| 10⁴ | 0.272106 | 0.087121 | 0.258312 | 0.013793 | 0.9569 | 0.0408 | 0.9372 | 0.086 | 0.003739 | 0.2725 |
| ∞ (BO limit) | 0.249338 | – | 0.249338 | 0 | 0.95865 | 0.03926 | 0.93966 | 0 | 0 | – |

The table uses the CODATA-2018 masses of the scan grid; the CODATA-2012 shift is far below the digits shown. All 37 masses, with errors, are in `table.csv`.

### 1. The R|ρ entanglement is not zero, and it does not vanish at M → ∞

The entanglement is of order 1 nat at every M. **It tends to a finite limit S(∞) = S_orient(∞) = 0.249338 nats.** This limit is rigorous: it is computed from the exact clamped-nuclei 1sσg function at R_e = 1.99719, expanded in partial waves about the midpoint, and converged to 1e−8 (w₀ = 0.95865, w₂ = 0.03926, w₄ = 1.84e−3, …).

The reason: in the BO limit Ψ → χ(R) φ(ρ; R⃗). The electron is in a pure state in the body frame, but in space-fixed coordinates its partial-wave content (l = 2, 4, … about the X–X axis) is locked to the direction R̂. So the remaining entanglement is purely **orientational (frame) entanglement**. It is not a breakdown of BO; it is the statement that only the total L = 0 is conserved, not the nuclear rotational angular momentum.

The radial part vanishes, as BO requires (§3).

### 2. Orientational entanglement and rotational decoherence

Trace out the electron. The X–X axis is then left with orientational coherence D(Θ) = Σ w_l P_l(cos Θ):
* D(90°) = 0.73 (M = 0.05), with a minimum of 0.7047 at M ≈ 0.475;
* D(90°) = 0.934 (H₂⁺), approaching 0.9397 as M → ∞.

Three points follow.
1. **D(180°) = 1 exactly at every M.** There are only even l, so the electron cannot distinguish R̂ from −R̂; this is the homonuclear exchange symmetry.
2. **The decoherence of the axis by its own electron saturates as M grows** (D → 0.94 rather than → 0). It is set by the anisotropy of σg at R_e, not by the mass.
3. **The internal electron does not localize the molecular axis.** Making the axis classical in the J = 0 ground state would require an external environment (collisions, radiation, a field). The "emergence of orientation" is therefore not internal. Conversely, the orientational entanglement is largest (and D smallest) near M ≈ 0.4–0.5, where the system is least molecule-like.

### 3. Large-M (BO) behaviour

**Fitted** (`fits.txt`; window M ≥ 500; local slopes over M ≥ 20):

| quantity | fitted exponent in κ | local slopes (approaching) | expected (derivation below) |
|---|---|---|---|
| S_orient − S_orient(∞) | 2.027 ± 0.003 | 2.13 → 2.015 | κ² |
| I(\|R\|;\|ρ\|) | 1.995 ± 0.000 | 2.11 → 1.996 | κ² |
| p₂ (second Schmidt coefficient, l = 0) | 1.957 ± 0.005 | 1.69 → 1.975 | κ² |
| E_N (radial) | 1.039 ± 0.001 | 1.11 → 1.034 | κ |
| S_rad | 1.68 (effective) | 1.56 → 1.71 | κ² ln(1/κ) |
| S − S(∞) | 1.83 (effective) | | mixture of κ² and κ² ln(1/κ) |

Fixed form for S_rad (M ≥ 206.77): S_rad = κ² (0.391 ln(1/κ) + 0.479), with residuals ≤ 1e−5.

**Derivation of the leading behaviour.** This is a justified leading-order argument, not a theorem. Write φ(ρ; R) ≈ φ(ρ; R_e) + (R − R_e) ∂_Rφ. The second Schmidt coefficient of each block is then p₂ ≈ ⟨δR²⟩ ‖∂_Rφ⊥‖², which scales as ΔR² ∝ κ². Consequences:
* the von Neumann entropy of a block with p₂ ≪ 1 is p₂(1 − ln p₂), i.e. κ² ln(1/κ);
* the log-negativity is ≈ 2√p₂, i.e. κ;
* the mutual information of a nearly Gaussian joint distribution is ≈ ½ corr², i.e. κ²;
* S_orient − S_orient(∞) comes from the R-dependence of w_l(R) averaged over δR, so it scales as ⟨R⟩ − R_e and ⟨δR²⟩, both ∝ κ².

All four predictions are consistent with the measured local slopes. The coefficient 0.391 of the κ² ln(1/κ) term is fitted; it is not the naive single-block value, because several l blocks contribute.

**BO product state against exact.**
* For M ≲ 20 the BO product state overestimates S strongly: 2.42 against 1.17 at M = 3, and 0.758 against 0.645 at M = 20. This is because the BO P(R) is much too broad there (cf. the Lindemann ratio in mass_scans).
* At large M the difference S_BO − S is 2.1e−3 (M = 1836), 1.0e−3 (4000) and 4.0e−4 (10⁴), roughly ∝ 1/M = κ⁴, i.e. adiabatic-correction order. This is from three points only.
* Part of this difference may come from linear interpolation of the DVR P_BO(R): the BO norm deviates from 1 by ≤ 5e−4, although each spectrum is normalized. So the κ⁴ attribution is tentative.

### 4. Small M and features near the markers

* **Small M.** S varies little for M ≤ 1 (1.495 at M = 0.05, 1.516 at Ps⁻). The M → 0 limit is the R|ρ entanglement of an H⁻-like ion with a fixed centre (finite); the trend at 0.05–0.1 shows it is approached slowly. This limit was not computed separately.
* **Maxima.** Every measure has one smooth maximum (minimum for w₀ and D) at small M. Locations come from a cubic spline in ln M through the scan grid plus 16 dense points on 0.25–0.6:

| measure | extremum | location M |
|---|---|---|
| E_N | max | ≈ 0.35 (jittery, see Validation) |
| I(\|R\|;\|ρ\|) | max | 0.357 |
| S_rad | max | 0.369 |
| S₀ | max | 0.387 |
| S | max (1.60403) | 0.414 |
| S_orient | max | 0.433 |
| S_lin | max | 0.445 |
| w₀ | min | 0.460 |
| D(90°) | min | 0.475 |

  The dense points are fixed-basis re-solves; at M = 0.4 and 0.5 they reproduce the scan values exactly. Locations are spline-resolved to about ±0.01. The quantities are flat (S changes by 3e−4 over 0.35–0.5), so the exact positions are sensitive to the descriptor and only modestly to the numerics.
* **Comparison with the markers.**
  * The entanglement maxima form a band at **0.35–0.48**. It lies between the Lindemann maximum (0.2–0.3) and King's g-minimum (0.4924) / the E/μ maximum (≈ 0.5), and below the CM-density crossover (0.635).
  * The radial measures peak lowest; the orientational ones peak highest, close to 0.49.
  * Nothing is sharp: no kink, no derivative jump, smooth spline curvature.
  * Inflection points at M ≈ 0.13–0.19 and 2.3–3.8 are smooth and grid-resolved.
  * There is nothing special at 0.543, 0.635, 1.11, ~1.4, 4.21, 9.53 or 11.49. The threshold masses (4.21, 9.53, 11.49) concern excited/rotational states and leave no trace in the ground-state entanglement, as expected for a ground state that stays bound and analytic in M.

  Conclusion: entanglement is another smooth, descriptor-dependent crossover marker near M ≈ 0.4–0.5, not a transition signal.

### 5. Literature

I know of **no published value** for the R|ρ (Jacobi) electron–nuclear entanglement of H₂⁺ or Ps⁻ that I can quote reliably. Papers on electron–nuclear or electron–electron entanglement in H₂⁺, H⁻ and Ps⁻ exist, but they mostly use different bipartitions (particle-based, or model Hamiltonians). From memory I cannot give numbers or references with confidence (unverified), so nothing is compared. A targeted literature search is needed before the paper claims novelty.

## Rigorous versus fitted

* **Rigorous** (exact identities, or exact results up to the stated quadrature/basis errors):
  * the block-diagonal Schmidt structure and the even-l-only property;
  * S = S_orient + S_rad, and D(180°) = 1;
  * every tabulated entropy (± 1e−5);
  * S(∞) = S_orient(∞) = 0.249338, w_l(∞) and D(90°, ∞) = 0.93966 (the infinite-mass limit of the exact wave function, assuming Ψ → BO product in the body frame, which is standard);
  * E_N is a lower bound on what survives the local angular traces.
* **Leading-order argument (justified, not proven):** exponents κ² (p₂, MI, S_orient − S∞), κ (E_N), κ² ln(1/κ) (S_rad).
* **Fitted:**
  * all quoted exponents and prefactors, and the coefficients 0.391 and 0.479;
  * extremum locations (spline);
  * S_BO − S ∝ 1/M (three points).

## Files

| file | content |
|---|---|
| `entcore.py` | engine: grids, partial-wave projection, spectra and split, radial negativity and mutual information, exact / N−4 / BO-product amplitudes, electronic M → ∞ weights |
| `validate.py`, `validation.txt` | tests above |
| `run.py` | resumable scan (modes exact, fine, prevN, bo, dense) |
| `results_{exact,fine,prevN,bo,dense}.csv` | raw rows per mass (S, S_lin, S_orient, S_rad, S₀, w₀…w₆, w_lmax, λ_max, K = 1/purity, E_N, negativity, MI, norms, timing) |
| `data/ent_<mode>_M*.npz` | l list, w_l, S_l, top Schmidt coefficients overall and per block, grid description |
| `analyze.py`, `table.csv`, `fits.txt` | errors (quadrature, basis, l-tail), BO columns, p₂, D(90°), fits, extrema |
| `fig1_entropies_vs_M.png` | S, S_orient, S_rad, S_lin, E_N, MI against M (exact solid, BO product dashed, markers in grey); w₀, w₂, D(90°), p₂ |
| `fig2_largeM_scaling.png` | log–log against κ, with κ and κ² guides |
| `fig3_orientation.png` | D(Θ) for selected M; w_l against l |
| `fig4_convergence.png` | quadrature, basis and l-tail errors |
| `log_*.txt` | run logs |

## Limitations

* Only the R|ρ split is computed. Other splits (for example particle-based) are kinematics-dominated and not comparable; see the harmonic_model README §5.
* E_N relies on a truncated local basis and is reliable only to about 5e−3 for M < 1.
* There is no independent check of the radial M → 0 limit.
* The BO-product comparison carries a P(R) interpolation error of about 5e−4.

## Body-frame check against Sanz-Vicario et al. PRA 96, 022503 (2017) (added 2026-11-13)
`bodyframe.py`: reduced density of the scalar |R| with the orientation traced out, Gamma = sum_l (2l+1) A_l A_l^T / tot. H2+ (M = 1836.15267343): S_vN = 0.047954 bits = 0.033239 nats, S_L = 0.010466, lambda_1 = 0.99474; base vs fine grid agree to 2e-8 bits. Sanz-Vicario et al. (their m_p = 1836.152673, energy 6.6e-5 above exact): 0.0479102 bits, S_L = 0.0104556. Agreement 0.1%. Results in `results_bodyframe.csv`.
