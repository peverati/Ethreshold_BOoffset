"""Validation of the unsymmetrized three-body solver for unequal nuclear masses: HD+ (J = 0, v = 0, 1).

Written 2026-11-16.  No crossing search; physical masses only (allowed before the prediction freeze).
Run from the project root (~5-15 min, < 3 GB):
    THREEBODY_BLAS_THREADS=8 python3 paper_data/unequal_masses/validate_hdplus.py
Output: validation_hdplus.txt.

References (both read from project files, not from memory):
  * Bubin & Adamowicz program paper (J. Chem. Phys. 152, 204102 (2020); masses from its example input inout_1_1500.txt:
    m_d = 3670.4829652, m_p = 1836.15267261): 1500 CECG, v=0 -0.597897968588, v=1 -0.589181829513
    (variational upper bounds).
  * Baskerville et al. 2016 (CODATA 2014 m_p = 1836.15267389,
    m_d = 3670.48296785): p+d+e- -0.59789796860 ("all digits converged").
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "8")
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
from threebody.coupled import CoupledBasis, block  # noqa: E402
from threebody.rotational import optimize_scales  # noqa: E402
from threebody.system import ThreeBodySystem  # noqa: E402

SETS = {"Bubin (ATOM-MOL-nonBO example masses)": ((1836.15267261, 3670.4829652), (-0.597897968588, -0.589181829513)),
        "Baskerville 2016 (CODATA 2014)": ((1836.15267389, 3670.48296785), (-0.59789796860, None))}
W = (1, 1, 2)
NS = (24, 28, 32, 36)


def main():
    out = []
    P = lambda s: (print(s, flush=True), out.append(s))
    P(f"# validate_hdplus.py {time.strftime('%Y-%m-%d %H:%M')}  particles (p, d, e) = (0, 1, 2), L = 0, no symmetry, weights {W}")
    (mp, md), _ = SETS["Bubin (ATOM-MOL-nonBO example masses)"]
    s = ThreeBodySystem((mp, md, 1.0), (1.0, 1.0, -1.0), name="HD+")
    t0 = time.perf_counter()
    al, e0, n = optimize_scales(s, 0, 20, None, state=0, x0=(9.36, 9.36, 1.41), weights=W)
    P(f"scales (N=20, root 0): {al}  E={e0:.12f}  ({n} evals, {time.perf_counter() - t0:.0f}s)")
    for N in NS:
        t1 = time.perf_counter()
        b = CoupledBasis(0, block(N, al, W), None)
        for name, ((mp, md), ref) in SETS.items():
            sys_ = ThreeBodySystem((mp, md, 1.0), (1.0, 1.0, -1.0))
            e, _ = b.solve(sys_, 2, vectors=False)
            d0 = e[0] - ref[0]
            d1 = (e[1] - ref[1]) if ref[1] is not None else float("nan")
            P(f"N={N:2d} size={b.size:5d}  {name:40s} E0={e[0]:.13f} (E0-ref {d0:+.2e})  E1={e[1]:.13f} (E1-ref {d1:+.2e})")
        P(f"   ({time.perf_counter() - t1:.0f}s; build {b.build_time:.0f}s)")
        del b
    P(f"E_th (p + D(1s)) = {-0.5 * md / (md + 1):.12f}")
    open(os.path.join(HERE, "validation_hdplus.txt"), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
