"""Critical masses for the lowest natural-parity L = 1 ('A') and L = 2 ('S') states of {X+, X+, e-}.

Threshold for every natural-parity L and both exchange symmetries: X + (Xe)(1s),
E_th = -M / (2 (M + 1)), reached with relative angular momentum lambda = L (README, section 10).

Method (fixed basis => mass-independent matrices, Brent root of E_0(M) - E_th(M)):
  * scales of a single compact block optimized for root 0 near the crossing;
  * optional diffuse blocks (scales (alpha_d, alpha_atom, alpha_atom) and mirror) for the
    X - (Xe) relative motion, alpha_atom = mu = M/(M+1) (exp(-mu (s1+s2)/2) = exp(-mu r12), the 1s atom X1 e);
  * crossing mass M*(basis) for a nested sequence of bases. Each M* is a rigorous upper bound
    to the exact critical mass (variational root below E_th proves binding), up to rounding;
  * binding energies E_th - E_0 at masses above M* for a linear threshold-law fit.

Usage: python examples/rotational_critical.py L sym [outdir]
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

from threebody import symmetric_family, family_threshold
from threebody.coupled import CoupledBasis, IdxBlock, block
from threebody.rotational import optimize_scales
from threebody.threshold import diffuse_index

L = int(sys.argv[1])
SYM = sys.argv[2]
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                         "paper_data", "rotational")
M_BO = {(1, "A"): 5.5374983484, (2, "S"): 12.7397228380}[(L, SYM)]
LOG = open(os.path.join(OUT, f"critical_L{L}{SYM}.txt"), "a")


def log(s):
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


def crossing(b, lo, hi, tol=1e-12):
    f = lambda M: b.solve(symmetric_family(M), 1, tol, vectors=False)[0][0] - family_threshold(M)
    flo, fhi = f(lo), f(hi)
    if not (flo > 0 > fhi):
        raise ValueError(f"no bracket: f({lo})={flo:.3e} f({hi})={fhi:.3e}")
    return brentq(f, lo, hi, xtol=1e-11, rtol=1e-13)


def diffuse_blocks(alpha_atom, alpha_d, n_c, n_d):
    idx = diffuse_index(n_c, n_d)[:, [1, 0, 2]]  # diffuse in s0: X0 far from the (X1 e) atom
    B = IdxBlock((alpha_d, alpha_atom, alpha_atom), idx, name=f"D{n_c}x{n_d}@{alpha_d:g}")
    return [B, B.mirror()]


def main():
    t_start = time.time()
    log(f"\n==== L={L} {SYM}  start {time.ctime()}  BO M_c = {M_BO}")
    # 1. scales: optimize root 0 at a mass where the state is clearly bound, then near crossing
    M_opt = M_BO
    al = None
    for it in range(3):
        al, e, n = optimize_scales(symmetric_family(M_opt), L, 12, SYM, x0=al)
        b = CoupledBasis(L, block(14, al), SYM)
        try:
            Mx = crossing(b, 0.5 * M_BO, 2.0 * M_BO)
        except ValueError as exc:
            log(f"  it{it}: {exc}")
            M_opt *= 1.3
            continue
        log(f"  it{it}: alphas opt at M={M_opt:.6f}: {al}; crossing (N=14) M*={Mx:.8f}")
        M_opt = Mx * 1.01
    res = dict(L=L, sym=SYM, alphas=al, compact={}, diffuse={})
    jpath = os.path.join(OUT, f"critical_L{L}{SYM}.json")
    if os.path.exists(jpath):
        old = json.load(open(jpath))
        if old.get("alphas") and np.allclose(old["alphas"], al):
            res["compact"] = {int(k): v for k, v in old["compact"].items()}
            log("  reusing compact crossings from json")
    # 2. nested compact sequence
    lo, hi = 0.5 * M_BO, 2.0 * M_BO
    for N in ((12, 14, 16, 18, 20, 22) if L == 1 else (12, 14, 16, 18, 20)):
        if N in res["compact"]:
            hi = res["compact"][N] * 1.0000001
            continue
        t0 = time.time()
        b = CoupledBasis(L, block(N, al), SYM)
        Mx = crossing(b, lo, hi)
        res["compact"][N] = Mx
        hi = Mx * 1.0000001
        log(f"  compact N={N} size={b.size} M*={Mx:.10f} ({time.time() - t0:.0f}s)")
        json.dump(res, open(os.path.join(OUT, f"critical_L{L}{SYM}.json"), "w"), indent=1)
    # 3. compact N=20 + diffuse blocks of increasing size
    mu = M_opt / (M_opt + 1)
    NB = 20 if L == 1 else 18
    hi = res["compact"][NB] * 1.0000001  # compact NB + diffuse contains compact NB
    blocks = [block(NB, al)]
    dlist = ((4, 12, 0.4), (4, 16, 0.15), (5, 16, 0.05), (5, 20, 0.02), (6, 20, 0.008))
    for (n_c, n_d, ad) in (dlist if L == 1 else dlist[:3]):
        t0 = time.time()
        blocks = blocks + diffuse_blocks(mu, ad, n_c, n_d)  # nested sequence
        b = CoupledBasis(L, blocks, SYM)
        tag = "+".join(B.name for B in blocks[1::2])
        try:
            Mx = crossing(b, lo, hi)
        except ValueError as exc:
            log(f"  diffuse {tag}: {exc}")
            continue
        pr = b.project(1e-12)
        res["diffuse"][tag] = Mx
        log(f"  compact{NB} + {tag}: size={b.size} kept={pr['kept']} M*={Mx:.10f} ({time.time() - t0:.0f}s)")
        hi = Mx * 1.0000001
        json.dump(res, open(os.path.join(OUT, f"critical_L{L}{SYM}.json"), "w"), indent=1)
    # 4. binding energies above M* with the largest basis (for a threshold-law fit)
    rows = []
    Mstar = min(list(res["compact"].values()) + list(res["diffuse"].values()))
    for x in (1.002, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2):
        M = Mstar * x
        e = b.solve(symmetric_family(M), 1, 1e-12, vectors=False)[0][0]
        rows.append((M, float(family_threshold(M) - e)))
        log(f"  M={M:.8f}  E0={e:.13f}  E_th-E0={family_threshold(M) - e:.4e}")
    res["binding"] = rows
    Ms, Eb = np.array(rows).T
    sel = Ms <= Mstar * 1.051
    p = np.polyfit(Ms[sel], Eb[sel], 2)
    roots = np.roots(p)
    log(f"  quadratic fit of binding (M <= 1.05 M*): coeffs {p}; roots {roots}")
    res["total_time_s"] = time.time() - t_start
    json.dump(res, open(os.path.join(OUT, f"critical_L{L}{SYM}.json"), "w"), indent=1)
    log(f"  done in {time.time() - t_start:.0f}s")


if __name__ == "__main__":
    main()
