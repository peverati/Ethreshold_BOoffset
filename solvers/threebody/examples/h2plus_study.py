"""H2+ (M = 1836.15267343): convergence with total degree N and anisotropic truncations.

Usage: python examples/h2plus_study.py      (writes results/h2plus.txt; several minutes)
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from threebody import symmetric_family, optimize_alphas, convergence_table, M_PROTON  # noqa: E402
from threebody.convergence import format_table  # noqa: E402

s = symmetric_family(M_PROTON)
ref = -0.597139063
out = []
for weights, Ns, Nopt in [((1, 1, 1), [16, 20, 24, 28, 32], 14), ((1, 1, 2), [24, 28, 32, 36, 40], 20)]:
    t0 = time.perf_counter()
    al, e, nev = optimize_alphas(s, Nopt, weights=weights)
    tab = convergence_table(s, al, Ns, nstates=4, weights=weights)
    r = tab["last"]
    txt = (f"=== H2+ weights={weights} (w0 a0 + w1 a1 + w2 a2 <= N), scales optimized at N={Nopt}\n"
           + format_table(tab, ref=ref)
           + f"\nroots below threshold at largest N: {r.n_below_threshold()}; time {time.perf_counter()-t0:.0f} s\n")
    print(txt, flush=True)
    out.append(txt)
with open(os.path.join(os.path.dirname(__file__), "..", "results", "h2plus.txt"), "w") as fh:
    fh.write("\n".join(out))
