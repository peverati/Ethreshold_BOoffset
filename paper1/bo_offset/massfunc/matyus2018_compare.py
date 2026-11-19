"""Compare our m_v(R), m_r(R) with Matyus, J. Chem. Phys. 149, 194111 (2018) (doi:10.1063/1.5050401).

Her printed numbers (no table of the functions in the paper; values are in her supplementary material, not available here):
  * Eq. (67): Cartesian body-fixed mass-correction tensor at rho = 3 bohr (z-axis embedding, nuclei 1, 2):
      same-nucleus  xx = yy = A = 0.322, zz = B = 0.462; cross-nucleus xx = yy = a = 0.178, zz = b = 0.038.
  * Eqs. (71), (84), (86): delta m_vib = B - b = 0.424, delta m_rot = A - a = 0.144 at rho = 3 bohr,
    entering as m_n -> m_n + delta m (per nucleus), i.e. identical to our M -> M + m.
  * Fig. 1 (read by eye, about +-0.01): delta m_vib maximum about 0.66 near 5-6 bohr; both curves level at 0.50.
Convention: delta m_vib = m_v, delta m_rot = m_r (ours). Schwenke EPAPS for the same R uses m = -C M^2 (Section 12).
Reads results_fine.csv and literature_data/schwenke2001_EPAPS_h2p_tab.txt; prints a short report. Light, no new calculation.
"""
import os
import numpy as np
import pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(H, "results_fine.csv")).sort_values("R")

# Schwenke EPAPS (same parsing as schwenke_epaps_compare.py)
rows = []
for line in open(os.path.join(H, "literature_data", "schwenke2001_EPAPS_h2p_tab.txt")).read().splitlines()[5:40]:
    p = line.split()
    if len(p) == 5:
        rows.append([float(x) for x in p])
t = np.array(rows)
M = 1836.15267
sR, s_mv, s_mr = t[:, 0], -t[:, 1] * 1e-7 * M**2, -t[:, 4] * 1e-9 * M**2 * t[:, 0] ** 2

R0 = 3.0
mv = float(np.interp(R0, d.R, d.m_v)); mr = float(np.interp(R0, d.R, d.m_r))
i = int(np.argmin(abs(sR - R0)))
A, B, a, b = 0.322, 0.462, 0.178, 0.038
print("Matyus 2018 vs ours at rho = R = 3 bohr")
print(f"  sum rules (translation carries the electron): A+a = {A+a:.3f}, B+b = {B+b:.3f} (exact 1/2)")
print(f"  {'':10s} {'Matyus':>8s} {'ours':>9s} {'diff':>8s} {'Schwenke':>9s}")
print(f"  {'m_v':10s} {B-b:8.3f} {mv:9.5f} {B-b-mv:+8.4f} {s_mv[i]:9.5f}   (Schwenke at R = {sR[i]})")
print(f"  {'m_r':10s} {A-a:8.3f} {mr:9.5f} {A-a-mr:+8.4f} {s_mr[i]:9.5f}")
print(f"  Matyus prints 3 decimals: rounding +-0.0005")
k = int(np.argmax(d.m_v.values))
print(f"ours: m_v maximum {d.m_v.values[k]:.4f} at R = {d.R.values[k]:.2f} bohr (Matyus Fig. 1: about 0.66 near 5-6 bohr)")
for Rx in (10.0, 15.0, 20.0, 40.0):
    print(f"ours at R = {Rx:5.1f}: m_v {np.interp(Rx, d.R, d.m_v):.4f}, m_r {np.interp(Rx, d.R, d.m_r):.4f}")
# Moss misprint cross-check: her Table I 'deviation' for (v,J) = (0,2) and (10,4) vs ours (Section 7)
print("Moss Table 4 misprints, deviation Moss - computed (cm^-1):")
print(f"  (0,2):  Matyus +0.059, ours {0.177-0.116:+.3f}")
print(f"  (10,4): Matyus -0.273, ours {1.037-1.308:+.3f}")
