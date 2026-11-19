"""Validation of the L = 3, 4 extension (angular.tensor_components) on H2+ rotational levels.

Written 2026-11-15.  M = 1836.15267247 (CODATA 2006), weights (1, 1, 2), fixed scales from the L = 1
optimization of examples/rotational_validation.py (9.95811007, 9.95811007, 1.42227562).

* L = 3 ('A'): (v=0, J=3) vs Tang, Yan, Shi, Babb, PRA 90, 012524 (2014), arXiv:1407.1973 Table 1
  (Tang et al., PRA 90, 012524 (2014), Table I): -0.5955576389803092(7).
* L = 4 ('S'): no verified reference in hand.  Check: E(0,J) for J = 0..3 (same table) fixed by
  E = E0 + B x - D x^2 + H x^3, x = J(J+1) (exact interpolation), extrapolated to J = 4. The neglected
  next term contributes L' x(x-2)(x-6)(x-12) = 40320 L' at x = 20; with L' ~ H^2/D ~ 6.5e-14 that is
  ~3e-9 Eh, so this check is good to a few 1e-9 only (result 2026-11-15: E(N=28) - interpolation = -3.1e-9,
  i.e. below it, which a variational energy can only be if the interpolation is off by at least that).

Usage: python examples/rotational_validation_L34.py [outdir]   (~10-20 min, < 3 GB)
"""
import os
import sys
import time

import numpy as np

from threebody import symmetric_family
from threebody.rotational import convergence_table_L, format_table_L

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                         "paper1", "thresholds_Lv")
M_P06 = 1836.15267247
AL = (9.95811007, 9.95811007, 1.42227562)
REF = {0: -0.597139063079392297758, 1: -0.5968737387847130778, 2: -0.5963452054891147, 3: -0.5955576389803092}


def main():
    x = np.array([J * (J + 1) for J in range(4)], float)
    c = np.linalg.solve(np.vander(x, 4, increasing=True), np.array([REF[J] for J in range(4)]))
    e4 = float(np.polyval(c[::-1], 20.0))
    s = symmetric_family(M_P06)
    with open(os.path.join(OUT, "validation_h2p_L34.txt"), "w") as fh:
        P = lambda t: (print(t, flush=True), fh.write(t + "\n"), fh.flush())
        P(f"# H2+ L=3,4 validation  {time.strftime('%Y-%m-%d %H:%M')}; M={M_P06}, weights (1,1,2), alphas {AL}")
        P(f"# J=4 reference from the J=0..3 Dunham-type interpolation: {e4:.13f} (coefficients E0,B,-D,H = {c})")
        for L, sym, ref in ((3, "A", REF[3]), (4, "S", e4)):
            t0 = time.time()
            tab = convergence_table_L(s, L, sym, AL, [20, 24, 28], nstates=1, weights=(1, 1, 2))
            P(format_table_L(tab, [ref]))
            P(f"L={L}: E(N=28) - ref = {tab['rows'][-1]['energies'][0] - ref:+.2e}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
