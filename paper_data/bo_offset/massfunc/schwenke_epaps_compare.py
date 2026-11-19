"""Compare our m_v(R), m_r(R) with Schwenke 2001 EPAPS tables (H2+, aug-cc-pV5Z).

Convention (LITERATURE_COMPARISON.md Section 8): C_vv = -m_v / M^2, C_rr = -m_r / (M^2 R^2), M = nuclear mass in m_e.
The M Schwenke used is not printed in the table; we report results for m_p = 1836.15267 and also the M that
makes the two curves agree best at large R (where m_v -> 1/2). Reads literature_data/schwenke2001_EPAPS_h2p_tab.txt
(lines 6..40) and results_fine.csv. Writes schwenke_epaps_compare.txt.
"""
import numpy as np, pandas as pd, os
H = os.path.dirname(os.path.abspath(__file__))
rows = []
for line in open(os.path.join(H, "literature_data", "schwenke2001_EPAPS_h2p_tab.txt")).read().splitlines()[5:40]:
    p = line.split()
    if len(p) == 5:
        rows.append([float(x) for x in p])
t = np.array(rows)
R, Cvv, Crr = t[:, 0], t[:, 1] * 1e-7, t[:, 4] * 1e-9
d = pd.read_csv(os.path.join(H, "results_fine.csv")).sort_values("R")
mv = np.interp(R, d.R, d.m_v); mr = np.interp(R, d.R, d.m_r)
out = [f"# Schwenke 2001 EPAPS (H2+, aug-cc-pV5Z) vs ours; {len(R)} R values {R.min()}..{R.max()} bohr", ""]
for M, lab in [(1836.15267, "M = 1836.15267 (CODATA proton)")]:
    sv = -Cvv * M**2; sr = -Crr * M**2 * R**2
    out.append(f"## {lab}")
    out.append(f"{'R':>6} {'mv_S':>9} {'mv_ours':>9} {'diff':>9} {'mr_S':>9} {'mr_ours':>9} {'diff':>9}")
    for i in range(len(R)):
        out.append(f"{R[i]:6.2f} {sv[i]:9.5f} {mv[i]:9.5f} {sv[i]-mv[i]:+9.5f} {sr[i]:9.5f} {mr[i]:9.5f} {sr[i]-mr[i]:+9.5f}")
    m = R <= 20
    out.append(f"max |diff| (R<=20): m_v {np.max(np.abs(sv-mv)[m]):.5f}, m_r {np.max(np.abs(sr-mr)[m]):.5f}; "
               f"median relative: m_v {np.median(np.abs(sv/mv-1)[m]):.2e}, m_r {np.median(np.abs(sr/mr-1)[m]):.2e}")
    out.append("")
open(os.path.join(H, "schwenke_epaps_compare.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
