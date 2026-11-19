"""Validation of the L = 1, 2 extension against verified literature values.

He (infinite nuclear mass): 2^1P (L=1, 'S'), 2^3P (L=1, 'A').
  Drake & Bondy, arXiv:2606.00814, Table 1 (checked online 2026-11-13):
  2^1P = -2.123843086498101359241(22), 2^3P = -2.133164190779283205140(8).
H2+ (m_p = 1836.15267247, CODATA 2006): (v, L) levels, Zhang, Yan et al., arXiv:1407.1973,
  Table 1 (checked online 2026-11-13; agree with Korobov):
  (0,1) -0.5968737387847130778(1), (1,1) -0.58690432091919159(5)   [L=1, 'A']
  (0,2) -0.5963452054891147(2),    (1,2) -0.5864036315280(8)       [L=2, 'S']
  (0,0) -0.597139063079392297758(4), (1,0) -0.58715567916469513(2)  [L=0, 'S']

Usage: python examples/rotational_validation.py he|h2p|all  [outdir]
"""
import json
import os
import sys
import time

import numpy as np

from threebody import helium_like, symmetric_family
from threebody.rotational import convergence_table_L, format_table_L, optimize_scales

OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                         "paper1", "rotational")
REF_HE = {"S": [-2.123843086498101], "A": [-2.133164190779283]}
M_P06 = 1836.15267247
REF_H2P = {0: [-0.597139063079392, -0.587155679164695], 1: [-0.596873738784713, -0.586904320919192],
           2: [-0.596345205489115, -0.586403631528]}


def log(fh, s):
    print(s, flush=True)
    fh.write(s + "\n")
    fh.flush()


def run_he():
    he = helium_like(2.0)
    with open(os.path.join(OUT, "validation_he.txt"), "w") as fh:
        log(fh, "He (m_nuc = inf), natural parity L = 1; scales optimized for root 0 at N = 12")
        for sym, name in (("S", "2^1P"), ("A", "2^3P")):
            t0 = time.time()
            al, e, n = optimize_scales(he, 1, 12, sym, x0=(1.0, 1.0, 2.0))
            log(fh, f"\n{name}: optimized alphas {al} E(N=12)={e:.13f} ({n} evals, {time.time() - t0:.0f}s)")
            tab = convergence_table_L(he, 1, sym, al, [12, 14, 16, 18, 20, 22], nstates=2, observables=True)
            log(fh, format_table_L(tab, REF_HE[sym]))
            log(fh, "observables (last N): " + json.dumps(tab["observables"][0]))


def run_h2p(Ls=(0, 1, 2), Ns=(24, 28, 32, 36), N_opt=20):
    s = symmetric_family(M_P06)
    with open(os.path.join(OUT, f"validation_h2p_L{''.join(map(str, Ls))}.txt"), "w") as fh:
        log(fh, f"H2+ M = {M_P06} (CODATA 2006), weights (1,1,2); scales optimized for root 0 at N = {N_opt}")
        for L in Ls:
            sym = "S" if L % 2 == 0 else "A"
            t0 = time.time()
            al, e, n = optimize_scales(s, L, N_opt, sym, x0=(9.36, 9.36, 1.41), weights=(1, 1, 2))
            log(fh, f"\nL={L} ({sym}): alphas {al} E(N={N_opt})={e:.13f} ({n} evals, {time.time() - t0:.0f}s)")
            tab = convergence_table_L(s, L, sym, al, list(Ns), nstates=2, weights=(1, 1, 2), observables=True)
            log(fh, format_table_L(tab, REF_H2P[L]))
            log(fh, "observables root0: " + json.dumps(tab["observables"][0]))


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("he", "all"):
        run_he()
    if what in ("h2p", "all"):
        run_h2p()
