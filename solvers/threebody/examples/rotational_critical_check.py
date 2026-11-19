"""Sensitivity of the critical-mass upper bound to the compact block size, with the diffuse blocks fixed.

Reads the scales from critical_L{L}{sym}.json (written by rotational_critical.py) and computes the
crossing mass M* for compact N in Ns plus the first n_diff diffuse blocks of that script.
Usage: python examples/rotational_critical_check.py L sym n_diff N1 N2 ... [--out dir]
"""
import json
import os
import sys
import time

sys.argv, OUTARG = (sys.argv[:sys.argv.index("--out")], sys.argv[sys.argv.index("--out") + 1]) \
    if "--out" in sys.argv else (sys.argv, None)
L, SYM, ND = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
NS = [int(x) for x in sys.argv[4:]]
sys.argv = [sys.argv[0], str(L), SYM] + ([OUTARG] if OUTARG else [])
import rotational_critical as rc  # noqa: E402  (opens its log in append mode)
from threebody.coupled import CoupledBasis, block  # noqa: E402

res = json.load(open(os.path.join(rc.OUT, f"critical_L{L}{SYM}.json")))
al = tuple(res["alphas"])
dlist = ((4, 12, 0.4), (4, 16, 0.15), (5, 16, 0.05), (5, 20, 0.02), (6, 20, 0.008))[:ND]
# same atom scale as rotational_critical.py: mu at M_opt = 1.01 x (crossing at N = 14 with the final scales)
M_opt = res["compact"]["14"] * 1.01
mu = M_opt / (M_opt + 1)
out = res.setdefault("compact_check", {})
for N in NS:
    t0 = time.time()
    blocks = [block(N, al)]
    for (n_c, n_d, ad) in dlist:
        blocks += rc.diffuse_blocks(mu, ad, n_c, n_d)
    b = CoupledBasis(L, blocks, SYM)
    Mx = rc.crossing(b, 0.5 * rc.M_BO, 1.5 * rc.M_BO)
    out[f"N{N}+D{ND}"] = Mx
    rc.log(f"  check: compact{N} + {ND} diffuse sets: size={b.size} kept={b.project(1e-12)['kept']} "
           f"M*={Mx:.10f} ({time.time() - t0:.0f}s)")
    json.dump(res, open(os.path.join(rc.OUT, f"critical_L{L}{SYM}.json"), "w"), indent=1)
