# threebody: deterministic variational solver for three-particle Coulomb S states

*Last updated: 2026-11-17 20:42*

> **Mass constants (2026-11-17):** Paper 1 quotes physical-mass values at CODATA 2022 (m_p/m_e = 1836.152673426, m_μ/m_e = 206.7682827); they are in `paper_data/codata2022/` (README, results.csv) and `solvers/threebody/results/codata2022/`. Examples and tables in this README that name other masses (CODATA 2018 1836.15267343, 206.768283; CODATA 2010 1836.15267245, 206.7682843, used for the King et al. comparison and called "CODATA 2012" in older notes) are kept as documentation of those runs; they differ from the CODATA 2022 values only in the last digits.

Non-relativistic, non-Born-Oppenheimer, total L = 0 states of three particles
with arbitrary masses and charges (in particular the {X+, X+, e-} family:
masses (M, M, 1), charges (+1, +1, -1)). Basis: Laguerre functions of the
perimetric coordinates (Pekeris-type), i.e. an orthogonal-polynomial rewriting
of a Hylleraas basis. All matrix elements are exact. There is no random
sampling, so for a fixed basis the results are a smooth, reproducible function
of M.

Python 3 with numpy, scipy and sympy (sympy only builds the operator tables;
tests also use it). Tested with the Anaconda Python at `python3`
(numpy 1.26.4, scipy 1.13.1, sympy 1.13.2). Nothing was installed. mpmath is
available but **not needed**: in double precision the conditioning is fine
(see §1).

```
threebody/
  laguerre.py     exact 1D integrals (Jacobi and derivative matrices), evaluation, mixed scales
  operators.py    S, T_i, V_ij, r_ij^n as polynomial term tables (sympy derivation)
  system.py       ThreeBodySystem, symmetric_family(M), helium_like(Z), thresholds
  basis.py        PerimetricBasis (truncation, exchange symmetry), matrix assembly
  solver.py       MatrixSet, solve(), Result (observables, g(r), save), overlap_between()
  convergence.py  optimize_alphas(), convergence_table(), extrapolate()
  scan.py         scan_masses(), CSV and g(r) writers
  __main__.py     command line interface
tests/            pytest suite, 38 tests (24 L = 0 core, 14 L = 1-4 in test_angular.py)
examples/         scripts that regenerate results/*.txt and *.csv
results/          validation.txt, h2plus.txt, scan_ground.csv, scan_excited1.csv
```

## 1. Formulation (and why)

### Coordinates
Particles 0, 1, 2 with masses m_i (infinity allowed) and charges q_i. For total
L = 0 the wave function depends only on the three distances r01, r02, r12. We use
**perimetric coordinates** (one per vertex of the particle triangle)

    s_k = r_kj + r_kl - r_jl ,   {k,j,l} = {0,1,2},   0 <= s_k < inf (independently),
    r_ij = (s_i + s_j) / 2,

    dtau = 8 pi^2 r01 r02 r12 dr01 dr02 dr12 = 8 pi^2 (r01 r02 r12 / 4) ds0 ds1 ds2 .

(8 pi^2 is dropped everywhere; it cancels in every ratio.) The triangle
inequalities become simple positivity constraints, so the integration domain is a
product of three half-lines.

### Basis (Pekeris-type, no random parameters)

    f_n(s; alpha) = L_n(alpha s) exp(-alpha s / 2)        (L_n = Laguerre polynomial)
    F_a = f_a0(s0; alpha0) f_a1(s1; alpha1) f_a2(s2; alpha2),   w0 a0 + w1 a1 + w2 a2 <= N

Default w = (1,1,1): total degree <= N. This spans **exactly** the Hylleraas space
{r01^i r02^j r12^k, i+j+k <= N} x exp(-c01 r01 - c02 r02 - c12 r12) with
alpha0 = c01 + c02, alpha1 = c01 + c12, alpha2 = c02 + c12, but the Laguerre
form is far better conditioned than monomials (scaled overlap condition number
~2e6 at N = 28, no basis function ever discarded).

Exchange of identical particles 0 and 1 is s0 <-> s1 (needs alpha0 = alpha1):
symmetry-adapted functions n_I (F_a +/- F_a') with a' = (a1, a0, a2) give
separate symmetric ('S') and antisymmetric ('A') blocks. Spin is not included:
'S'/'A' is the spatial symmetry (for two electrons, S = singlet, A = triplet; for
two identical spin-1/2 X, the same).

There are only two (symmetric) or three (general) nonlinear parameters alpha.
They are chosen **deterministically** by Nelder-Mead minimization of the chosen
root at a modest size N_opt (default 12), started from a hydrogenic rule
(c_ij = |Z_ij| mu_ij for attractive pairs). The convergence table then uses these
fixed alphas at increasing N, which gives a strictly nested sequence of spaces and
hence strictly monotonic variational energies (Hylleraas-Undheim-MacDonald).

**Why this formulation**
* every matrix element is *exact* (no quadrature, no random numbers): the result is
  a smooth, reproducible function of the masses for fixed alphas;
* the Hamiltonian matrix is linear in the inverse masses and the couplings,
  H = sum_i (1/m_i) T_i + sum_ij Z_ij V_ij, with mass-independent T_i, V_ij: a mass
  scan at fixed basis only needs re-diagonalization, and dH/dM is available exactly;
* a single integer N controls the basis, nested spaces give clean convergence
  tables;
* the same machinery evaluates <r_ij^n> (n >= -1) exactly and overlaps between
  wave functions computed with different alphas (different M) in fixed
  physical coordinates.

### Hamiltonian for arbitrary masses (S states)
For f = f(r01, r02, r12) and particle i with partners j, k:

    grad_i f = (df/dr_ij) rhat_ij + (df/dr_ik) rhat_ik,
    cos_i = rhat_ij . rhat_ik = (r_ij^2 + r_ik^2 - r_jk^2) / (2 r_ij r_ik),
    Lap_i f = f_,ij,ij + (2/r_ij) f_,ij + f_,ik,ik + (2/r_ik) f_,ik + 2 cos_i f_,ij,ik .

H = -sum_i Lap_i / (2 m_i) + sum_{i<j} Z_ij / r_ij  (Z_ij = q_i q_j by default; any
pair can be switched off). Since f depends only on distances it is translation
invariant, so the centre-of-mass kinetic energy is zero automatically. Matrix
elements use the symmetric (integrated-by-parts) form

    <a|T|b> = sum_i 1/(2 m_i) int [a_ij b_ij + a_ik b_ik + cos_i (a_ij b_ik + a_ik b_ij)] dtau ,

where cos_i r01 r02 r12 = (r_ij^2 + r_ik^2 - r_jk^2) r_jk / 2 is a polynomial; with
d/dr_ij = d/ds_i + d/ds_j - d/ds_k every integrand (overlap, kinetic, Coulomb,
r_ij^n with n >= -1) becomes a polynomial in s times products of 1D basis
functions or their first derivatives.

**Verification** (tests/test_hamiltonian.py): (i) the distance-form Laplacian is
compared with a brute-force sympy Cartesian Laplacian in 9 coordinates for a
generic non-symmetric test function at random configurations (agreement 1e-11
relative) for each particle; (ii) the engine's S, V_ij and T_i matrices (first-
derivative form) are compared with direct Gauss-Laguerre quadrature of
<F_a| -Lap_i F_b /2> (second-order form), <F_a|F_b>, <F_a|1/r_ij|F_b> for a
non-symmetric basis with three different alphas (agreement ~1e-10 relative);
(iii) exact limits: e-e interaction switched off with an infinitely heavy nucleus
gives E = -Z^2 exactly (also with unequal light masses: -Z^2 (m0+m1)/2), the
mass-scaling identity E(lambda m) = lambda E(m) holds to 1e-11, and the symmetry-
adapted code agrees with the general unsymmetrized code after relabelling the
particles of an unequal-mass system (M, M, 1) to 1e-11.

### Exact 1D integrals
All integrals reduce to products of
M[a,b] = int (d^da f_a)(d^db f_b) s^p ds, da, db in {0,1}, obtained from the integer
Jacobi matrix of the Laguerre polynomials (int f_a f_b s^p ds = alpha^-(p+1) (J^p)_ab)
and the half-integer derivative matrix (d f_n/ds = alpha sum_k D_kn f_k,
D_kn = -1 for k<n, -1/2 for k=n). No cancellation, no extended precision needed.
The generalized eigenproblem is solved with canonical orthogonalization (S scaled
to unit diagonal; directions with eigenvalue < 1e-13 of the largest are dropped:
none were dropped in any run reported here).

## 2. How to run

```bash
cd solvers/threebody
PY=python3
$PY -m pytest -q tests                                   # full test suite, about 20 s
$PY -m threebody solve --M 1 --N 20 --nstates 4          # Ps-, scales optimized at N=12
$PY -m threebody converge --M 1836.15267245 --N 24 28 32 36 40 --weights 1 1 2 --N-opt 20 --nstates 4
$PY -m threebody solve --masses 1 1 inf --charges -1 -1 2 --identical --N 20     # He
$PY -m threebody scan --M 1 2 5 10 20 --N 16 20 24 --nstates 4 --out scan.csv \
      --save-dir scan_out --grid ../../templates/grid_0.05_to_25.txt          # CSV + npz + g(r)
$PY -m threebody scan --M 5 10 --N 16 20 24 --opt-state 1 --nstates 3 --out exc.csv   # scales tuned for root 1
```
Options: `--symmetry S|A|none|auto`, `--alphas a0 a1 a2` (fixed scales, no
optimization), `--N-opt` (size used for the scale optimization), `--weights`
(anisotropic truncation w0 a0 + w1 a1 + w2 a2 <= N). The BLAS thread count is
limited to 4 (override with `THREEBODY_BLAS_THREADS`), because OpenBLAS with all
10 cores was about 10x *slower* here for these matrix sizes.

Python API:
```python
from threebody import *
s = symmetric_family(1836.15267245)
al, _, _ = optimize_alphas(s, 20, weights=(1, 1, 2))
tab = convergence_table(s, al, [28, 32, 36, 40], nstates=4, weights=(1, 1, 2))
res = tab["last"]                       # Result at the largest N
res.energies, res.threshold, res.n_below_threshold()
res.observables(0)                      # <r_ij>, <r_ij^2>, <1/r_ij>, <1/r_ij^2>, T, V, virial, Lindemann
res.g_of_r(0, 1, r_grid)                # pair distribution, int 4 pi r^2 g dr = 1
res.dE_dinverse_mass(0)                 # Hellmann-Feynman <T_0> = dE/d(1/m_0) at fixed basis
res.save("H2p.npz")                     # basis definition + coefficients + energies (+ JSON metadata)
overlap_between(res_M, res_Mprime)      # <psi(M)|psi(M')> in fixed physical coordinates
```
For mass scans, `MatrixSet(basis)` can be built once and passed to
`solve(system, basis, matrices=ms)` for every M (only the diagonalization is
repeated). The saved `.npz` contains `tensors[state][a0,a1,a2]` (coefficients on
the primitive product functions), `alphas`, `N`, `weights` and the conventions as a
JSON string, so wave functions can be rebuilt or overlapped later without this
code. `overlap_between` handles different alphas and sizes exactly: the mixed-scale
1D integrals are polynomials times exponentials, done with Gauss-Laguerre rules.
Eigenvector signs are fixed by the convention <F_000|psi> > 0 (or by the largest
component when that overlap is below 1e-6). For excited states near avoided
crossings, use `overlap_between` to align states across M.

## 3. Observables
* `<r_ij^n>` for n = -1, 1, 2 (and any integer n >= -1): exact, by contracting the
  coefficient tensor with 1D moment matrices.
* `<1/r_ij^2>`: from the pair distribution by Gauss-Laguerre quadrature in r_ij
  (exact for alpha_i = alpha_j: the integrand is then a polynomial times
  exp(-(alpha_i+alpha_j) r)). The same routine reproduces the exact `<r>` and
  `<1/r>` to 1e-12 (tested).
* T, V, E = T + V and the virial ratio -V/T. Because alpha is optimized and a
  common rescaling of all alphas is a coordinate dilation, -V/T = 2 is satisfied to
  the optimizer tolerance: 1e-9 for light systems, 4e-8 for H2+ at default
  settings. This checks the optimization and the operators, not the convergence.
* Lindemann ratio of the X-X distance, sqrt(<r01^2> - <r01>^2) / <r01>.
* Pair distribution P_ij(R) (int P dR = 1) and g_ij(R) = P/(4 pi R^2): fix
  r_ij = R (s_i + s_j = 2R). The integral over the free perimetric coordinate is exact;
  the one over the finite segment uses Gauss-Legendre.
* `dE_dinverse_mass(i)` = <T_i>: the exact derivative of the variational eigenvalue
  with respect to 1/m_i at a fixed basis (Hellmann-Feynman; tested against finite
  differences). For the family, dE/dM = -(<T_0> + <T_1>)/M^2.

## 4. Convergence control and error estimates
`N` (total Laguerre degree) is the single size parameter. `convergence_table`
keeps the alphas fixed, so the spaces are nested and every root decreases
monotonically. `extrapolate` assumes a geometric tail over successive table entries: d_n = E_{n-1} - E_n,
rho = d_n/d_{n-1}, E_ext = E_n - d_n rho/(1-rho). The error is max(tail, change of
E_ext from the previous step). It is a heuristic estimate, not a bound. Flags:
`rho>=1` (not yet in the convergent regime), `nonmonotonic` (changes at the rounding
level, see below), `saturated` (change below ~1e-14).
**Rounding floor:** the scaled overlap has a condition number of about 2e6 at N = 28,
so eigenvalues carry rounding errors of about 1e-13 to 1e-11. Changes below that
level, and the `nonmonotonic` flags at M = 50 and 100 for state 1 (|change| < 2e-11),
are noise and not convergence information.

**Counting bound states.** Each variational root is an upper bound to the
corresponding exact eigenvalue (Hylleraas-Undheim-MacDonald). A root below E_th
therefore *proves* (up to rounding, margin 1e-10) that at least that many bound
states of that symmetry exist: `n_below_threshold()`. This is only a lower bound.
It is also capped by `nstates` (the CSV column `n_below_is_capped_by_nstates` flags
this). A root above threshold proves nothing: it may be discretized continuum, or a
weakly bound state that the basis does not yet describe. Its extrapolated value is
not meaningful when the sequence is accumulating at E_th.

## 5. Validation (results/validation.txt, results/h2plus.txt)

Targets were quoted from memory by the requester. Default total-degree basis.
Scales were optimized for the ground state at N = 12 and then kept fixed. "size" is
the number of symmetric basis functions.

| system | N (size) | E (variational) | extrapolated ± est. | target | E − target |
|---|---|---|---|---|---|
| He, m_nuc = ∞ | 28 (2360) | −2.90372437702875 | −2.90372437703255 ± 4e−12 | −2.903724377034 | +5.3e−12 |
| H⁻, m_nuc = ∞ | 28 (2360) | −0.52775101646188 | −0.52775101652007 ± 6e−11 | −0.527751016544 | +8.2e−11 |
| Ps⁻ (M = 1) | 28 (2360) | −0.26200507020127 | −0.26200507022359 ± 2e−11 | −0.262005070233 | +3.2e−11 |
| H₂⁺, M = 1836.15267245, total degree | 32 (3417) | −0.59713906306119 | −0.5971390630866 ± 3e−11 | −0.597139063 | (agrees to all 9 quoted digits) |
| H₂⁺, weights (1,1,2) | 40 (3311) | −0.59713906307925 | −0.59713906307925 ± 3e−13 | −0.597139063 | (agrees to all 9 quoted digits) |

In every case the variational value lies above the target and the extrapolated
value agrees with the target within its error estimate. (For H₂⁺ the target has
only 9 digits.) I recall a longer literature value, −0.597139063123, but that is
from memory and **unverified**. If it is correct, our value, internally converged
to about 1e−12, is 4e−11 too high. A different proton mass would not explain this:
dE/dM ≈ 1e−6 per electron mass, so the shift would be negligible. **This is an
open discrepancy at the 1e−10 level.** Check it against a trusted reference
before relying on digits beyond the tenth.
**Mass constants (2026-11-13):** every physical-mass value in this README (H₂⁺, Mu₂⁺) uses CODATA 2012, m_p = 1836.15267245 and m_μ = 206.7682843 (as King et al. 2013); the Tang et al. validation (§ L > 0) keeps their own m_p = 1836.15267247. Values were recomputed in `paper_data/codata2012/`; fixed-basis sequences are shifted by the first-order mass effect (Mu₂⁺ −6.04e−11, H₂⁺ +1.62e−12). The saved files in `results/` (`h2plus.txt`, `validation.txt`, `scan_*.csv`) still hold the earlier CODATA-2018 values. For H₂⁺ L=1,2 only 10 digits are quoted because the recomputed values come from an extrapolation of error up to 4.5e−11.
**Resolved (2026-11-13):** the remembered value is
correct for the older proton mass 1836.152701 (this solver gives −0.5971390631232 there,
`paper_data/king2013`). The argument above was wrong: the mass difference of 2.8e−5 m_e times
dE/dM ≈ 1.5e−6 gives exactly the 4e−11 shift. See also the L > 0 validation section below.
Other internal checks (all pass): exact separable limits; mass-scaling identity;
relabelling invariance; Cartesian-Laplacian and quadrature verification of the
operators; monotonic convergence; virial ratio; Hellmann-Feynman versus finite
differences; normalization and moments of g(r); He expectation values versus
well-known values (<r1> = 0.929472295, <r12> = 1.422070255, <1/r12> = 0.945818448,
<1/r1²> = 6.0174089, also quoted from memory; agreement 1e−8 or better).
Known excited levels come out right: He 2¹S −2.145974 (3¹S converging toward
−2.0613); H₂⁺ v = 1 −0.587155679. These are also from memory, not independently checked.

### Convergence tables (ground state, scales fixed)

He (alphas 2.8486, 2.8486, 4.8006)
```
 N  size   E0                   dE0
12   252  -2.90372437521442
16   525  -2.90372437678321   1.57e-09
20   946  -2.90372437697749   1.94e-10
24  1547  -2.90372437701808   4.06e-11
28  2360  -2.90372437702875   1.07e-11
```
Ps- (alphas 0.39587, 0.39587, 0.78418)
```
12   252  -0.26200505940958
16   525  -0.26200506896995   9.56e-09
20   946  -0.26200506993393   9.64e-10
24  1547  -0.26200507014319   2.09e-10
28  2360  -0.26200507020127   5.81e-11
```
Mu2+ (M = 206.7682843, CODATA 2012; alphas 4.9204, 4.9204, 1.4059). There is no literature
target; the value is internally converged.
```
12   252  -0.58512601860249
16   525  -0.58512609774690   7.91e-08
20   946  -0.58512609825192   5.05e-10
24  1547  -0.58512609825255   6.33e-13
28  2360  -0.58512609825255   (rounding level)
```
H2+ (M = 1836.15267245), total degree (alphas 8.0296, 8.0296, 1.2530, optimized at N = 14) and
anisotropic truncation a0 + a1 + 2 a2 <= N (alphas 9.3616, 9.3616, 1.4108, optimized at N = 20):
```
total degree                               weights (1,1,2)
 N  size  E0                 dE0           N  size  E0                 dE0
16   525  -0.59713244875271               24   819  -0.59713905530895
20   946  -0.59713869065246  6.2e-06      28  1240  -0.59713906291141  7.6e-09
24  1547  -0.59713904797472  3.6e-07      32  1785  -0.59713906307409  1.6e-10
28  2360  -0.59713906246705  1.4e-08      36  2470  -0.59713906307887  4.8e-12
32  3417  -0.59713906306119  5.9e-10      40  3311  -0.59713906307925  3.9e-13
```
H2+ lowest four S levels (weights (1,1,2), N = 40): −0.59713906307925,
−0.58715567909559, −0.5777519044, −0.5689084987 (roots 2 and 3 not recomputed at the CODATA-2012 mass, so quoted to 10 digits; their shifts are of order 1e−12). Estimated errors are
3e−13, 2e−12, 2e−11 and 6e−10. All four lie far below E_th = −0.49972784, so H2+
has at least 4 bound S states. Only 4 roots were requested; more roots would show
more vibrational levels.

## 6. Accuracy versus M (results/scan_ground.csv)
Default settings: total degree, N = 16, 20, 24, 28, scales optimized for the ground
state at N = 12 (warm-started along the scan). Each mass takes about 25 s, most of
it the scale optimization.

| M | E0 (N=28) | est. error E0 | E_th | E1 − E_th | roots < E_th (of 4) |
|---|---|---|---|---|---|
| 0.5 | −0.174185832909 | 3e−11 | −0.1666667 | +2.8e−4 | 1 |
| 1 (Ps⁻) | −0.262005070201 | 2e−11 | −0.25 | +4.0e−4 | 1 |
| 2 | −0.352687353126 | 9e−12 | −0.3333333 | +5.1e−4 | 1 |
| 5 | −0.450780241114 | 1e−12 | −0.4166667 | +5.9e−4 | 1 |
| 10 | −0.501800153396 | 1e−14* | −0.4545455 | +6.5e−5 | 1 |
| 20 | −0.536021015599 | 1e−13* | −0.4761905 | −8.7e−3 | 2 |
| 50 | −0.563750580187 | 1e−13* | −0.4901961 | −3.0e−2 | 3 |
| 100 | −0.576442412864 | 1e−13* | −0.4950495 | −4.7e−2 | 4 (capped) |
| 206.768 (Mu₂⁺) | −0.585126098252 | 1e−13* | −0.4975935 | −6.1e−2 | 4 (capped) |
| 500 | −0.591760947317 | 1e−12 | −0.4990020 | −7.5e−2 | 4 (capped) |
| 1000 | −0.595093299672 | 7e−10 | −0.4995005 | −8.2e−2 | 4 (capped) |
| 1836.15 (H₂⁺) | −0.597139059343 | 4e−8 | −0.4997278 | −8.8e−2 | 4 (capped) |

\* At the rounding floor: the realistic accuracy there is about 1e−12.

**Claimed accuracy of the ground-state energy** (default total-degree basis, N <= 28):
better than 1e−10 for M <= ~5, about 1e−12 for 10 <= M <= 500, about 1e−9 at
M = 1000, and only about 4e−8 at M = 1836 (H₂⁺). With the anisotropic truncation
(weights 1,1,2; N = 36–40) H₂⁺ reaches about 1e−12 internally (§5). Use it for
M >~ 500. The evidence is the monotonic, geometrically converging variational
sequences, the agreement with the literature targets at M = 1, ∞ and 1836, and
virial ratios equal to 2 within 1e−9 (4e−8 for H₂⁺ at default settings). Excited
states converge more slowly. Deeply bound vibrational levels at M >= 50 are as
accurate as the ground state. Levels within ~1e−3 of threshold are not converged
(errors 1e−5 to 1e−3) and need scales tuned for them (`--opt-state`).

**Second bound S state (symmetric).** With ground-state scales the second root lies
below threshold for M = 20 (first seen there) but not for M <= 10. With scales
optimized for root 1 (results/scan_excited1.csv):

| M | E1 (N=24, scales for root 1) | E_th | E1 − E_th |
|---|---|---|---|
| 2 | −0.3332139399 | −0.3333333333 | +1.2e−4 |
| 5 | −0.4165768194 | −0.4166666667 | +9.0e−5 |
| 10 | −0.4545839321 | −0.4545454545 | **−3.8e−5** |
| 20 | −0.4848548570 | −0.4761904762 | −8.7e−3 |

So for M = 10 a second symmetric S bound state is **proven**: a variational upper
bound lies 3.8e−5 below threshold. Its energy is not converged; the upper bound
improved by 1e−4 from the ground-state-optimized scales. For M <= 5 (and Ps⁻) no
second bound state is found, which is consistent with the belief that Ps⁻ has a
single bound state but does **not prove** it. The critical mass where the second
state appears is therefore between ~5 (not detected) and 10 (proven). It was not
located more precisely. Doing so needs a dedicated threshold-oriented basis
(see §9).

## 7. Extension to L > 0 (design notes; now implemented, see §10)

What would be needed (natural-parity states, L = 1, 2, ...):
* Wave function in the vector-coupled solid-harmonic form (Schwartz; Drake)
  Psi_LM = sum_{l1 + l2 = L} Y^{LM}_{l1 l2}(r_02, r_12) f_{l1}(s0, s1, s2),
  i.e. L+1 coupled radial components, each expanded in the same
  `PerimetricBasis` (possibly different N or alphas per component).
* Matrix elements between components (l1, l2) and (l1', l2'): the angular
  integrals give polynomials in cos theta_{02,12} = (r02^2 + r12^2 - r01^2)/(2 r02 r12)
  and the prefactors r02^{l1} r12^{l2} add polynomial factors; the kinetic energy
  acquires first-derivative cross terms and centrifugal-like terms. After
  multiplying by the measure r01 r02 r12 all integrands are again
  polynomial(s) x (first derivatives of 1D functions), so the existing engine
  (`PerimetricBasis._raw_matrix`, which takes any list of terms
  (coef, (p0,p1,p2), k, l)) evaluates them exactly without change.
* Code changes: (1) operators.py: functions returning term lists for each pair of
  components (derive them with sympy exactly as `kinetic_terms` does now and
  verify against a Cartesian Laplacian of Y^{LM} f, like tests/test_hamiltonian.py);
  (2) a `CoupledBasis` holding one PerimetricBasis per component and assembling the
  block matrix; (3) exchange symmetry maps component (l1, l2) to (l2, l1) with a
  phase (-1)^{l1 + l2 - L}; (4) observables: the radial-distance observables keep the
  same contraction code applied per component pair (block-diagonal in components
  for operators that depend only on distances, after angular integration).
  Unnatural parity (e.g. 1,2-even P^e states) needs the r02 x r12 factor in addition.
* The rest (solver, convergence, scan, CLI) is L-agnostic: `solve` only needs
  matrices S and H, and `Result` only needs a way to turn coefficients into
  per-component coefficient tensors.

## 8. Limitations
* Only L = 0 in the original modules (natural-parity L = 1-4: §10 and §11; unnatural parity not implemented). No spin; 'S'/'A' are spatial exchange symmetries of particles 0 and 1.
* One set of exponential scales per coordinate. Very diffuse near-threshold states
  converge slowly. A second ("double") set of scales would cure this; it is not
  implemented.
* The large-M (Born-Oppenheimer-like) regime needs the anisotropic truncation or
  large N. Matrix size and build time grow as N³ and N⁶ respectively: about 25 s
  at size 3400, memory about 1 GB.
* The scale optimization (Nelder-Mead, xatol 1e−4 in log alpha) is deterministic
  but stops at a finite tolerance. With warm starts the optimized alphas, and hence
  E at finite N, depend smoothly on M but also on the scan path. **For derivatives
  in M, keep the basis fixed** (pass `alphas=` / `--alphas`, or reuse a `MatrixSet`).
  Then E(M) is analytic and dE/dM is given exactly by `dE_dinverse_mass`.
* The extrapolation error is a heuristic estimate, not a bound.
* `n_below_threshold` gives a lower bound on the number of bound states only.

## 9. What was verified and what was not
Verified (automated tests in tests/, or the example scripts):
- Hamiltonian/Laplacian versus Cartesian sympy; matrix elements versus brute-force quadrature.
- Exact limits: separable, mass scaling, relabelling.
- Literature ground states of He, H⁻, Ps⁻, H₂⁺ (targets from memory).
- Monotonicity and virial ratio.
- Normalization of g(r); <r>, <1/r> from g(r) versus exact; overlaps (self = 1,
  different-scale overlaps, orthogonality of different roots).
- Hellmann-Feynman versus finite difference; save round trip; CSV scan.

Not verified / open:
- Mu₂⁺ and other intermediate-M energies: no reference values were compared.
- The H₂⁺ 1e−10-level question in §5.
- Antisymmetric ('A') states: the code path is exercised only by the relabelling
  test, which checks that the lowest A root lies above the lowest S root and that
  S ∪ A equals the unsymmetrized calculation. No literature value was checked.
- Systems without identical particles (e.g. HD⁺): exercised only by tests, not
  validated against literature.
- Accuracy of `<1/r²>` for alpha_i != alpha_j: the quadrature is not exact there;
  it converges with the node count, but this was not studied.

## Suggestions for next steps
1. For large M: use weights (1,1,2), or a dedicated basis whose X-X factor is
   localized (for example Laguerre in (r01 − R0)/sigma, or a James-Coolidge/Kolos-
   Wolniewicz-type elliptic basis with a Gaussian-like R factor). Alternatively
   two or more scale sets.
2. For thresholds: optimize the scales on root k, add a diffuse second scale set,
   and track E_k − E_th versus N and M. Fit the approach to threshold to locate the
   critical M.
3. Resolve the H₂⁺ 4e−11 question against a trusted reference value.
4. ~~Implement L = 1, 2 following §7.~~ Done (natural parity), §10.


## 10. Natural-parity states with L = 1, 2 (added 2026-11-13)

New files only: `threebody/angular.py` (operator term lists), `threebody/coupled.py`
(`CoupledBasis`), `threebody/rotational.py` (scale optimization, convergence tables,
observables), `tests/test_angular.py` (8 tests), `examples/rotational_validation.py`,
`examples/rotational_scan.py`, `examples/rotational_critical.py`,
`examples/rotational_compare.py`. A generated cache `threebody/_cache/angular_terms_v1.pkl`
holds the sympy-derived term lists (safe to delete; rebuilt in ~5 s for L = 1, ~40 s
for L = 2). Existing modules are unchanged except this README; the default S-state
behaviour is untouched. Results and logs: `paper_data/rotational/`.

### Formulation
Relative vectors a = x0 − x2, b = x1 − x2 (particle 2 = electron for {X,X,e}, nucleus for He).
A natural-parity state (parity (−1)^L) is written in the vector-coupled bipolar-harmonic form
used for two-electron atoms (Bhatia & Temkin, Rev. Mod. Phys. 36, 1050 (1964); Drake,
Phys. Rev. A 18, 820 (1978)):

    Psi_LM = sum_{l1+l2=L} |a|^l1 |b|^l2 [Y_l1(a^) x Y_l2(b^)]_LM  f_{l1 l2}(r01, r02, r12)

L + 1 components; each radial function is expanded in the same perimetric Laguerre
products (one or more blocks of scales). The components are coded as Cartesian irreducible
tensors (L=1: a_k, b_k; L=2: symmetric traceless parts of a a, a b, b b). All operators are
scalars, so the matrix element summed over the Cartesian tensor index is rotation invariant;
it is reduced to the three distances by evaluating it in a canonical frame (sympy, exact
rational arithmetic), times the usual measure r01 r02 r12. Kinetic energy in the
symmetric first-derivative form, sum_K grad_i(A_K F) . grad_i(A_K G). Every integrand is
again polynomial(s) x (first derivatives of 1D Laguerre functions): the existing exact 1D
machinery (`threshold.raw_matrix`, `laguerre`) is reused unchanged. H = sum (1/m_i) T_i +
sum Z_ij V_ij with mass-independent matrices (a mass scan only re-diagonalizes).
The completeness of the l1 + l2 = L form is taken from the literature above (not
re-derived); the numerical evidence is the agreement with verified He 2P and H2+ (v, L = 1, 2)
energies to 1e-10 ... 1e-13 below.

### Exchange symmetry and thresholds (derived, not assumed)
* P01 (exchange of the identical particles 0, 1) maps a <-> b, s0 <-> s1 and component
  (l1, l2) -> (l2, l1) with phase (−1)^(l1+l2−L) = +1 for natural parity. Symmetry-adapted
  functions n (A_c F ± A_Pc PF); for a P-invariant operator the matrix element is
  2 n n' (<A_c F|O|A_d G> ± <A_c F|O|A_Pd PG>). Both P = ±1 exist for every L.
* Which P belongs to the rotational tower of the ground vibronic manifold: in the BO limit
  Psi ~ phi_el(r; R) chi(R) Y_LM(R^), R = x0 − x1, r = electron position from the X–X
  midpoint. P01 sends R -> −R and leaves r unchanged. For 1s sigma_g, phi_el(r; −R) =
  phi_el(r; R) (gerade: the electronic function is symmetric under swapping the two
  centres), and Y_LM(−R^) = (−1)^L Y_LM(R^). Hence **P01 = (−1)^L: L = 0, 2 symmetric ('S'),
  L = 1 antisymmetric ('A')**, total parity (−1)^L (natural). The other symmetry at the
  same L correlates with 2p sigma_u (P01 = −(−1)^L), a repulsive curve with only a shallow
  long-range well. With spin: for spin-1/2 X (protons) 'S' pairs with the nuclear-spin
  singlet (para), 'A' with the triplet (ortho); for spinless bosonic X the 'A' levels
  (odd L of the ground manifold) would not exist.
  Numerical confirmation: at M = 1836, L=1 'A' root 0 is the (v=0, L=1) level (−0.5968737),
  far below E_th, while L=1 'S' is not attempted (diffuse 2p sigma_u); at M <= 10 the L=1 'S'
  roots lie above E_th (scan_L.csv).
* Threshold: the lowest continuum of natural parity L (either P) is X + (Xe)(1s) with
  relative angular momentum lambda = L (parity (−1)^lambda), so E_th = −M/(2(M+1)) for all
  L, as for S states (`coupled_threshold`). Both P are allowed asymptotically (the far X can
  be either particle). Unnatural parity would need an excited atom (n = 2). For lambda >= 1
  the open channel has a centrifugal barrier and a −alpha/(2R^4) polarization tail; the
  zero-energy solution u ~ R^−lambda is normalizable, so near M_c the binding vanishes
  ~ linearly and the state does not become a halo (unlike L = 0).

### Usage
```python
from threebody import symmetric_family, helium_like
from threebody.coupled import CoupledBasis, block
from threebody.rotational import optimize_scales, convergence_table_L, format_table_L
s = symmetric_family(1836.15267247)
al, _, _ = optimize_scales(s, 1, 20, "A", x0=(9.4, 9.4, 1.4), weights=(1, 1, 2))
tab = convergence_table_L(s, 1, "A", al, [24, 28, 32, 36], nstates=2, weights=(1, 1, 2))
print(format_table_L(tab))
b = CoupledBasis(2, block(20, al), "S")       # mass-independent matrices
e, C = b.solve(symmetric_family(50.0), nstates=2)
```
`CoupledBasis(L, blocks, symmetry)` accepts several `IdxBlock`s with different scales (e.g.
diffuse blocks for near-threshold states, as in `examples/rotational_critical.py`).
`symmetry=None` gives the unsymmetrized basis for non-identical particles (all T_i, V_ij stored).

### Validation (paper_data/rotational/validation_*.txt)
References checked online on 2026-11-13. He: G. W. F. Drake and A. T. Bondy,
arXiv:2606.00814, Table 1 (infinite nuclear mass). H2+: arXiv:1407.1973 (Hylleraas,
fully non-adiabatic), Table 1, m_p = 1836.15267247 (CODATA 2006); it agrees with
Korobov's values listed there. Scales optimized for root 0 at a small N, then fixed (nested
spaces, monotonic convergence).

| system, state | L, P | basis (size) | E (variational) | E − ref | extrap. − ref |
|---|---|---|---|---|---|
| He 2¹P | 1, S | N=22 (2300) | −2.12384308639634 | +1.0e−10 | +2.0e−11 (±8e−11) |
| He 2³P | 1, A | N=22 (2300) | −2.13316419070850 | +7.1e−11 | +1.4e−11 (±6e−11) |
| He 3¹P (root 1) | 1, S | N=22 | −2.05514634990126 | +1.2e−8 | −1.2e−9 (±1.3e−8) |
| He 3³P (root 1) | 1, A | N=22 | −2.05808106620289 | +1.8e−8 | −1.6e−9 (±2e−8) |
| H2+ (v,L)=(0,0) | 0, S | w(1,1,2) N=36 (2470) | −0.59713906307899 | +4.0e−13 | +2.6e−13 |
| H2+ (1,0) | 0, S | same | −0.58715567909077 | +7.4e−11 | +7.1e−11 |
| H2+ (0,1) | 1, A | w(1,1,2) N=36 (4750) | −0.59687373878466 | +5.7e−14 | +3.1e−14 |
| H2+ (1,1) | 1, A | same | −0.58690432091814 | +1.0e−12 | +5.6e−13 |
| H2+ (0,2) | 2, S | w(1,1,2) N=32 (5202) | −0.59634520548881 | +3.0e−13 | −7e−14 |
| H2+ (1,2) | 2, S | same | −0.58640363151524 | +1.3e−12 | +1.9e−12 |

All variational values lie above the references, as they must. Virial ratio −V/T = 2 to
5e−10 (He) and 1e−11 (H2+). The H2+ (0,0) agreement (2.6e−13) also **resolves the open
question of §5**: the remembered −0.597139063123 is the correct value for an older proton
mass, m_p = 1836.152701 (this solver gives −0.5971390631232 there; `paper_data/king2013`); the
verified reference at the CODATA-2006 mass 1836.15267247 is −0.597139063079392. The mass
matters at this level: the 5e−10 relative difference between 1836.15267247 and the
1836.15267343 (CODATA 2018, used in an earlier version of §5) shifts E by about −1.5e−12 (dE/dlnM ≈ −0.0028). §5 and the tables below now use the CODATA-2012 mass 1836.15267245, at which this solver gives −0.5971390630793, within 1e−13 of the Tang et al. value.
For L = 2 the canonical orthogonalization drops 15–47 near-dependent directions (cutoff
1e−12) at N = 24–32 with weights (1,1,2); the results above show this is harmless.

Tests (`tests/test_angular.py`, 8 tests, ~3 min, mostly sympy): L=0 term lists identical
to `operators.py`; kinetic integrand versus a brute-force Cartesian gradient (9 coordinates,
sympy) at random configurations for every particle and component pair, L = 1, 2 (agreement
1e−10 relative); S, T_i, V_ij matrices versus Gauss–Laguerre quadrature of the
*second-order* form with a Cartesian Laplacian (checks the integration by parts) for L = 1, 2
with three different scales; S ∪ A spectrum = unsymmetrized spectrum; hermiticity and
linearity in 1/M.

### Rotational critical masses (paper_data/rotational/critical_*.txt, *.json)
Mass at which the lowest L = 1 ('A') and L = 2 ('S') state first lies below E_th = −M/(2(M+1)).
Method: fixed basis (mass-independent matrices), Brent root of E_0(M) − E_th(M). Basis: one
compact total-degree block (scales optimized for root 0 near the crossing) plus nested
diffuse blocks (scales (alpha_d, mu, mu) and mirrors, alpha_d = 0.4 ... 0.008; one X far
from the (X e) 1s atom). Each crossing M* is a **rigorous upper bound** to M_c (a variational
root below E_th proves a bound state; rounding ~1e−13 Eh ⇔ ~3e−11 in M). The central value
is a geometric extrapolation in the compact size at fixed diffuse blocks (not a bound).

| state | compact only, N = 12 … 22 | + diffuse sets | compact N + diffuse (fixed) | best proven bound | extrapolated M_c | BO | BO + adiabatic |
|---|---|---|---|---|---|---|---|
| L=1, A | 4.21319, 4.21114, 4.21013, 4.20955, 4.20920, 4.20897 | N=20 + 1…5 sets: 4.2091246, 4.2083324, 4.2082849, 4.2082831, 4.2082830 | 4 sets, N=16/18/20/22: 4.20828396, 4.20828340, 4.20828312, 4.20828297 | 4.2082830 | 4.2082827 ± 5e−7 | 5.5374983 | 4.5705653 |
| L=2, S | 11.49003, 11.48971, 11.48958, 11.48952, 11.48948 (N ≤ 20) | N=18 + 1…3 sets: 11.4894301, 11.4894221, 11.4894221 | 3 sets, N=16/18/20: 11.48942332, 11.48942208, 11.48942181 | 11.4894218 | 11.4894217 ± 3e−7 | 12.7397228 | 11.7786164 |

BO values: `bo_h2plus` node counting (lowest level of rotation J bound below its own BO
threshold; two step sizes, Richardson; error < 1e−10). The binding just above M_c is linear in
M − M_c (slopes 3.1e−3 and 2.5e−3 Eh per electron mass for L = 1, 2), as expected for a
centrifugal barrier (lambda = L >= 1); the quadratic fit printed at the end of
`rotational_critical.py` is not used (it is biased by the far points).
Both L = 1 and L = 2 of the ground manifold are bound only for M > 11.4894: below that the
"rotational tower" does not exist (L = 1 alone for 4.208 < M < 11.489; neither below 4.208,
which includes Ps⁻-like M = 1). For comparison (paper_data/RESULTS_SUMMARY.md) the first
vibrational excitation (L = 0, 2nd 'S' state) appears at M = 9.5256.

### Mass scan of the rotational levels (paper_data/rotational/scan_L.csv, rotational_table.txt, fig_rotational.png)
Lowest two roots of L=0 'S', L=1 'A', L=2 'S' (and L=1 'S' for M <= 10). Scales: optimized
for L=0 and L=1 root 0 (N_opt = 12, or 16 with weights (1,1,2) for M > 150); L=2 reuses the
L=1 scales. Sequences N = 14, 18, 21 (M < 150) or 24, 28, 32 with weights (1,1,2).
d1 = E(L=1) − E(L=0), d2 = E(L=2) − E(L=0) (hartree); BO and BO+ad from `bo_h2plus`.

| M | E(L=0) | E(L=1) | E(L=2) | d1 | d2/d1 | M·d1 | E(v=1)−E(v=0) | d1(BO+ad)/d1 − 1 |
|---|---|---|---|---|---|---|---|---|
| 5 | −0.450780241108 | −0.419889382466 | (unbound) | 3.0891e−2 | – | 0.1545 | (v=1 unbound) | −4.9e−3 |
| 7 | −0.477943605538 | −0.449447840196 | (unbound) | 2.8496e−2 | – | 0.1995 | (unbound) | +4.3e−3 |
| 10 | −0.501800153396 | −0.477287816888 | (unbound) | 2.4512e−2 | – | 0.2451 | (bound, see §6) | +5.0e−3 |
| 20 | −0.536021015599 | −0.519901362845 | −0.496129580099 | 1.6120e−2 | 2.4747 | 0.3224 | 5.117e−2 | +3.1e−3 |
| 50 | −0.563750580187 | −0.555852034437 | −0.542038371483 | 7.8986e−3 | 2.7489 | 0.3949 | 4.321e−2 | +1.3e−3 |
| 100 | −0.576442412864 | −0.572140917518 | −0.564121202390 | 4.3015e−3 | 2.8644 | 0.4302 | 3.467e−2 | +6.9e−4 |
| 206.768 | −0.585126098252 | −0.582928452634 | −0.578683577569 | 2.1977e−3 | 2.9316 | 0.4544 | 2.630e−2 | +3.4e−4 |
| 1836.153 | −0.597139063045 | −0.5968737388 | −0.5963452055 | 2.6532e−4 | 2.9920 | 0.4872 | 9.983e−3 | +4e−5 |

Estimated errors (extrapolation): bound levels 1e−15 … 1e−12 for 5 <= M <= 207, ~5e−11 at
M = 1836 (scan settings; the validation runs reach 1e−13). Levels marked unbound lie above
E_th in every basis used, and the M = 1–3 rows of scan_L.csv (all L > 0 roots above E_th) are
discretized continuum: their "spacings" mean nothing. For comparison, at M = 1836 BO gives
d2/d1 = 2.9920 as well, and M·d1 = 0.4872 equals 2<1/R²>(BO, v=0) = 2·29.156 cm⁻¹·M/219474.6
= 0.4878 to 0.1%.
Observations (not over-interpreted): (i) M·d1 rises monotonically from 0.15 (M = 5) to 0.49
(M = 1836): the spacing is ∝ 1/M only asymptotically; (ii) d2/d1 approaches the rigid-rotor 3
from below (2.47 at M = 20, 2.99 at 1836), and the BO value follows the exact one closely
(2.4638 vs 2.4747 at M = 20): the deviation from 3 is mostly centrifugal distortion, already
present in BO; (iii) the BO+adiabatic d1 deviates from the exact one by ≈ +0.07/M (relative)
for M >= 20, consistent with an O(1/M) non-adiabatic correction; plain BO crosses from +0.8%
(M = 20) to −1e−4 (M = 207). (iv) The tower {L = 0, 1, 2} exists only for M > 11.4894.


## 11. L = 3, 4 (added 2026-11-15)

`angular.tensor_components(c)` gives the rank-L symmetric traceless Cartesian tensor of component c for any L
(harmonic projection of prod_k (u_k . x), one entry per index multiset with multinomial weight: 10 entries for
L = 3, 15 for L = 4); the term lists use the weighted sums. L <= 2 is unchanged (same flat lists, cache valid).
The L = 3, 4 terms are in `_cache/angular_terms_v1.pkl` (one-time sympy cost ~15 min for L = 4).
Tests (6 new, all pass; full suite 38 passed in 9.5 min): weighted contraction = full 3^L contraction from an
independent detracer formula (L = 3, 4); S + A = unsymmetrized spectrum (L = 3, 4); pointwise kinetic integrand
(L = 3); second-order-form quadrature (L = 3).
Validation (`examples/rotational_validation_L34.py`, output `paper_data/thresholds_Lv/validation_h2p_L34.txt`):
H2+ (0, 3), N = 28: E - ref = +1.2e-11 (extrapolated +1.2e-12) vs Tang et al. 2014; (0, 4): -3.1e-9 below a
J = 0..3 Dunham interpolation of the same table, consistent with the neglected x^4 term (no verified literature
value for (0, 4) in hand). Basis sizes: L = 3 two channels, L = 4 three channels per block (S or A).
