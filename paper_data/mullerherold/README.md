# Cross-check of solvable model B against Muller-Herold (JCP 124, 014105, 2006)

*Last updated: 2026-11-13 19:36*

Files: `derive_mh.py` (sympy re-derivation), `compare.py` (table, quadrature, pointwise density check; imports
`solvers/harmonic_model` read-only). Run with `python3`.

## 1. Re-derivation of his critical relation (Sec. III)
Jacobi: r1 = q2-q1, r2 = q3-(q1+q2)/2, r3 = CM; q1 = -r1/2 - (m_u/M) r2, M = 2m_e+m_u, mu = 2 m_e m_u/M.
Ground state: |psi1|^2 = r1^{2t} exp(-sqrt(m_e) r1^2), t(t+1)=m_e (checked: both the r^{t-2} and r^{t+2}
terms of the radial equation vanish); |psi2|^2 = exp(-sqrt(8 mu) r2^2) (note: the *density* carries sqrt(8mu), not 8mu).
Density of a pair particle about the CM: n(x) = Int d3r1 r1^{2t} e^{-c r1^2} exp(-eps2 (x + r1/2)^2),
with c = sqrt(m_e) and eps2 = sqrt(8 mu)(M/m_u)^2 = 4 sqrt(m_e) (M/m_u)^{3/2}.
Angular average to O(x^2): n(x)/n(0)-1 ~ x^2 <-eps2 + eps2^2 r1^2/6>, weight r1^{2t+2} e^{-(c+eps2/4) r1^2}, so
`<r1^2> = (t+3/2)/(c+eps2/4)`. Curvature zero iff eps2 t/6 = c (his "t = 6 s", s = sqrt(t(t+1))/eps2; reproduced).
Inserting eps2: (M/m_u)^{3/2} = 3/(2t), i.e.

**Corrected Eq. 3.13:  m_u = 2 m_e / [ (3/2)^{2/3} (sqrt(m_e+1/4) - 1/2)^{-2/3} - 1 ]**   (t = sqrt(m_e+1/4)-1/2).

Checks: m_u = 200 -> m_e = 3.4573 (his 3.46, OK). Curvature positive (bimodal) for m_e above the line, negative below.
**Discrepancy:** as m_u -> inf, M/m_u -> 1, t -> 3/2, so the asymptote is m_e = t(t+1) = **15/4 = 3.75**
(m_e = 3.743 at m_u = 1e4). His text and Fig. 1 caption say 15/16; that cannot be an asymptote of his own relation
(m_e = 15/16 is crossed already at m_u = 2.17, and m_e = 3.46 at m_u = 200 exceeds it). Likely a typo/misprint (15/4); not verifiable
beyond the PDF text. Leave a flag in the paper.

## 2. Mapping to model B
MH = model B with k = K = 4 (k' = 0), g = 1, pair mass M = m_e, third mass m_u. Our code fixes third mass = 1; rescaling
lengths by 1 and masses by m_u (energies by 1/m_u) gives **M = m_e/m_u, g = m_u, k = K = 4 m_u** (wrapper `mapB` in compare.py;
package untouched). Algebraically our condition lam(gM) = (3/2) sqrt(K/k)(2M+1)^{-3/2} with lam(lam+1) = gM = m_e (so lam = t) and K/k = 1
becomes (M/m_u)^{3/2} = 3/(2t) with M/m_u -> (2M+1): identical to MH. Only g M and K/k matter, so the default
model-B call `critical_mass_cm(g=1)` (k=K=1) **is** the MH model at m_u = 1: M_c = m_e = 0.63306 (not "a different number", as the older note said; it is his
curve evaluated at m_u = 1). General statement: the critical-mass condition of model B is MH's Eq. 3.13 with K/k = 1; for K/k != 1 replace 3/2 by (3/2) sqrt(K/k).

## 3. Comparison (compare.py output)
| m_u | m_e (MH Eq. 3.13 corrected) | m_e (model B, critical_mass_cm via mapping) | m_e (independent 2D quadrature + FD curvature) |
|---|---|---|---|
| 0.5 | 0.4373223003 | 0.4373223003 | 0.437324 |
| 1 | 0.6330580666 | 0.6330580666 | 0.633060 |
| 5 | 1.3786880565 | 1.3786880565 | 1.378695 |
| 20 | 2.3060974115 | 2.3060974115 | 2.306114 |
| 200 | 3.4573248884 | 3.4573248884 | 3.457361 |
| 1836.15 | 3.7138578474 | 3.7138578474 | - |
| 1e4 | 3.7432712438 | 3.7432712438 | - |

MH vs model B: relative difference <= 1.6e-14 (machine precision). Quadrature (direct integration of the MH wavefunctions on the z axis,
Richardson-extrapolated second difference) agrees to ~1e-5 (finite-difference limited). Pointwise check: package `heavy_density_cm` under the mapping reproduces the
quadrature n(z)/n(0) to all printed digits (10) at (m_e=2, m_u=10, z=0.1,0.3,0.6). Sign of `cm_curvature_coefficient` flips across the line as expected.
Conclusion: model B at (k=K, any g) reproduces MH exactly; the two papers' CM-density transitions are the same criterion on the same model.
Note: the model-B bond length R_e = (4g/K)^{1/4} equals 1 here (his r1 at the classical minimum).
