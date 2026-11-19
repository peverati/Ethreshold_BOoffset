# CODATA-2012 copy of examples/validation_tables.py (generated 2026-11-14); writes into this folder.
"""Regenerate the validation and convergence tables quoted in README.md.

Usage:  python examples/validation_tables.py [Nmax]     (writes results/validation.txt)
Run time: a few minutes with Nmax = 28.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np  # noqa: E402

from threebody import helium_like, symmetric_family, optimize_alphas, convergence_table  # noqa: E402
M_PROTON, M_MUON = 1836.15267245, 206.7682843  # CODATA 2012 (Paper 1)
from threebody.convergence import format_table  # noqa: E402

Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 28
Ns = list(range(12, Nmax + 1, 4))
cases = [
    ("He (Z=2, m_nuc=inf)", helium_like(2.0), -2.903724377034),
    ("H- (m_nuc=inf)", helium_like(1.0), -0.527751016544),
    ("Ps- (M=1)", symmetric_family(1.0), -0.262005070233),
    ("M=10", symmetric_family(10.0), None),
    ("Mu2+ (M=206.7682843)", symmetric_family(M_MUON), None),
    ("H2+ (M=1836.15267245)", symmetric_family(M_PROTON), -0.597139063),
]
out = []
for name, s, ref in cases:
    t0 = time.perf_counter()
    al, e, nev = optimize_alphas(s, 12)
    tab = convergence_table(s, al, Ns, nstates=3)
    txt = f"=== {name}: scales optimized for the ground state at N=12 ({nev} evals)\n" + format_table(tab, ref=ref)
    r = tab["last"]
    obs = r.observables(0)
    txt += (f"\nN={Ns[-1]}: virial -V/T = {obs['virial_ratio']:.10f}; roots below threshold = {r.n_below_threshold()}"
            f"; S cond (scaled) = {r.info['S_cond_scaled']:.2e}; dropped = {r.info['dropped']}"
            f"\n<r01> = {obs['r01']:.10f}  <r02> = {obs['r02']:.10f}  Lindemann(01) = {obs['lindemann01']:.6f}"
            f"\ntotal time {time.perf_counter() - t0:.1f} s\n")
    print(txt, flush=True)
    out.append(txt)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "validation.txt"), "w") as fh:
    fh.write("\n".join(out))
