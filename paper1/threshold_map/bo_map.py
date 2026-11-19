"""BO and BO+ad threshold map of the 1s sigma_g manifold of {X+, X+, e-}: every (J, v) with M_c <= M_MAX.

Written 2026-11-16 (threshold_map, step 1).  Run from the project root (~5-10 min, cheap):
    python3 paper1/threshold_map/bo_map.py

Method: zero-energy coupling-constant problem of paper1/bo_offset/coupling.py (lam = M + 1/2, any L;
c = 0 pure BO, c = 1 BO + adiabatic correction counted at its own asymptote), as in
paper1/thresholds_Lv/bo_brackets.py.  The roots lam_c(L, v) in (0.5, M_MAX + 0.5] are located on a dense lam
grid (step 0.25) and refined with brentq (coupling.critical).
Cross-checks: (i) pure BO roots against bo_h2plus count_levels bisection (Numerov, ppw 16/32 Richardson) for a
selection of (J, v); (ii) BO and BO+ad level counts at M = 100 from count_levels against the map.
2p sigma_u: first BO threshold is J=0 at M = 330.6 (thresholds_Lv/bo_brackets.txt), so no 2p sigma_u level
enters below M ~ 100.

Exchange symmetry of a 1s sigma_g level: (-1)^J  (odd J 'A', even J 'S'); in the exact three-body problem
(J, v) is root k = v of the natural-parity L = J sector of that symmetry.

Output: bo_map.json, bo_map.txt (also the order of appearance and the comparison with Frolov 1992 II).
"""
import json
import math
import os
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("BO_H2PLUS_THREADS", "4")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "paper1", "bo_offset"))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
import numpy as np  # noqa: E402

import coupling as C  # noqa: E402

M_MAX = 105.0
LMAX = 10
# exact values known before this project (do not recompute)
EXACT = {(0, 1): 9.525634, (0, 2): 29.63939, (0, 3): 59.31866, (0, 4): 98.26550,
         (1, 0): 4.2082827, (1, 1): 17.094282, (2, 0): 11.4894217, (3, 0): 21.376199, (4, 0): 33.914349}
# Frolov, J. Phys. B 25, 3059 (1992) (paper II): quoted y_t and the smallest mass at which his table binds the level
FROLOV = {(0, 1): (10.045, 10.05), (0, 2): (30.70, 30.75), (0, 3): (61.135, 61.15), (1, 1): (17.125, 17.1),
          (1, 2): (39.375, 39.45), (2, 1): (27.715, 27.75), (2, 2): (53.295, 53.30)}
FROLOV_ORDER = ["(0,0)", "(1,0)", "(0,1)", "(2,0)", "(1,1)", "(2,1)", "(0,2)", "(1,2)", "(2,2)", "(0,3)", "(1,3)"]
# his list "1S, 3P, 2 1S, 1D, 2 3P, 2 1D, 3 1S, 3 3P, 3 1D, 4 1S, 4 3P" (no F, G): n^{2S+1}L -> (J=L, v=n-1)

OUT = []


def P(s=""):
    print(s, flush=True)
    OUT.append(s)


def roots(L, c):
    lo = 0.5 + 1e-6
    hi = M_MAX + 0.5
    n = int((hi - lo) / 0.25) + 1
    return [float(x) - 0.5 for x in C.critical(L, c, lo, hi, n=n)]


def bo_h2plus_crit(J, n, Mc):
    """Pure BO critical mass of level v = n-1 by count_levels bisection (ppw 16 and 32, Richardson)."""
    from bo_h2plus import count_levels

    def crit(ppw, lo, hi, rtol=1e-11):
        cnt = lambda M: count_levels(M, J, adiabatic=False, ppw=ppw)
        assert cnt(lo) < n <= cnt(hi), (J, n, lo, hi, cnt(lo), cnt(hi))
        while hi / lo - 1 > rtol:
            mid = math.sqrt(lo * hi)
            if cnt(mid) >= n:
                hi = mid
            else:
                lo = mid
        return math.sqrt(lo * hi)
    M1 = crit(16.0, Mc * 0.97, Mc * 1.03)
    M2 = crit(32.0, M1 * 0.99, M1 * 1.01)
    return M2 + (M2 - M1) / 15, abs(M2 - M1) / 15


def main():
    t0 = time.time()
    P(f"# BO / BO+ad threshold map, 1s sigma_g manifold of {{X+,X+,e-}}  ({time.strftime('%Y-%m-%d %H:%M')})")
    P(f"# M_c = lam_c - 1/2 (electron masses), all roots with M_c <= {M_MAX}; L = J = 0..{LMAX}")
    res = {}
    for L in range(LMAX + 1):
        for c, name in ((0, "BO"), (1, "BO+ad")):
            r = roots(L, c)
            res[f"L={L} {name}"] = r
            P(f"L={L:2d} {name:6s} M_c(v=0,1,..) = " + (", ".join(f"{x:.8f}" for x in r) or "none"))
    # ---- table of levels
    levels = []
    for L in range(LMAX + 1):
        bo = res[f"L={L} BO"]
        # BO+ad roots: the J=0 v=0 level is bound for every lam > 1/2 with the adiabatic correction, so its root is
        # absent; label each BO+ad root by v = number of BO roots below it (offsets ~1 << level spacing)
        ad = {sum(1 for b in bo if b < a): a for a in res[f"L={L} BO+ad"]}
        assert len(ad) == len(res[f"L={L} BO+ad"]), (L, bo, ad)
        for v in range(max(len(bo), max(ad) + 1 if ad else 0)):
            levels.append(dict(J=L, v=v, sym="A" if L % 2 else "S",
                               M_BO=bo[v] if v < len(bo) else None, M_BOad=ad.get(v),
                               exact=EXACT.get((L, v)), frolov=FROLOV.get((L, v), (None, None))[0],
                               frolov_bound_at=FROLOV.get((L, v), (None, None))[1]))
    # BO+ad of every level with BO <= M_MAX: if BO+ad is missing it lies above M_MAX (BO+ad < BO always so far)
    levels.sort(key=lambda d: d["M_BOad"] if d["M_BOad"] is not None else 1e9)
    P("\n# Levels ordered by BO+ad M_c (sector: exchange symmetry (-1)^J, root k = v)")
    P("| order | (J,v) | sector | M_c BO | M_c BO+ad | BO - BO+ad | exact (known) | Frolov 1992 II y_t (bound at) | status |")
    P("|---|---|---|---|---|---|---|---|---|")
    for i, d in enumerate(levels):
        st = "done" if d["exact"] else ("TARGET" if d["M_BOad"] is not None and d["M_BOad"] <= 101.0 else "out of range")
        if d["J"] == 0 and d["v"] == 0:
            st = "always bound"
        d["status"] = st
        f = lambda x, p=6: "" if x is None else f"{x:.{p}f}"
        diff = f(d["M_BO"] - d["M_BOad"], 4) if (d["M_BO"] is not None and d["M_BOad"] is not None) else ""
        fr = f"{d['frolov']} ({d['frolov_bound_at']})" if d["frolov"] else ""
        P(f"| {i + 1} | ({d['J']},{d['v']}) | L={d['J']} '{d['sym']}' root {d['v']} | {f(d['M_BO'])} | {f(d['M_BOad'])} | "
          f"{diff} | {f(d['exact'])} | {fr} | {st} |")
    targets = [d for d in levels if d["status"] == "TARGET"]
    P(f"\n# {len(targets)} targets (BO+ad <= 101, no exact value yet): "
      + ", ".join(f"({d['J']},{d['v']})" for d in targets))
    # ---- order of appearance vs Frolov
    order_ad = [f"({d['J']},{d['v']})" for d in levels if d["M_BOad"] is not None]
    order_bo = [f"({d['J']},{d['v']})" for d in sorted(levels, key=lambda d: d["M_BO"] if d["M_BO"] is not None else 1e9)
                if d["M_BO"] is not None]
    P("\n# Order of appearance (BO+ad): " + " < ".join(order_ad))
    P("# Order of appearance (BO):    " + " < ".join(order_bo))
    sub = [x for x in order_ad if x in FROLOV_ORDER]
    P("# Frolov 1992 II order (S,P,D only, his list): " + " < ".join(FROLOV_ORDER))
    P("# BO+ad order restricted to his levels:        " + " < ".join(sub)
      + ("   -> SAME ORDER" if sub == [x for x in FROLOV_ORDER if x != "(0,0)"] else "   -> DIFFERENT (see above)"))
    # ((0,0) has no BO+ad root: bound for every M; it is excluded from the comparison)
    # ---- cross-checks with bo_h2plus
    P("\n# Cross-check of pure BO roots with bo_h2plus count_levels bisection (ppw 16/32 Richardson)")
    chk = {}
    for (J, v) in ((1, 2), (2, 1), (5, 0), (6, 0), (7, 0), (2, 3), (1, 4)):
        bo = res[f"L={J} BO"]
        if v >= len(bo):
            continue
        val, err = bo_h2plus_crit(J, v + 1, bo[v])
        chk[f"({J},{v})"] = val
        P(f"  (J={J}, v={v}): bo_h2plus {val:.8f} (+/- {err:.1e});  coupling {bo[v]:.8f};  diff {val - bo[v]:+.1e}")
    from bo_h2plus import count_levels
    P("# Level counts at M = 100 (count_levels, ppw 16) vs map")
    cnt = {}
    for J in range(LMAX + 1):
        nb = count_levels(100.0, J, adiabatic=False, ppw=16.0)
        na = count_levels(100.0, J, adiabatic=True, ppw=16.0)
        mb = sum(1 for x in res[f"L={J} BO"] if x < 100.0)
        ma = sum(1 for x in res[f"L={J} BO+ad"] if x < 100.0) + (1 if J == 0 else 0)  # + always-bound (0,0)
        cnt[J] = (nb, na)
        P(f"  J={J}: count_levels BO {nb}, BO+ad {na};  map BO {mb}, BO+ad {ma}  {'OK' if (nb, na) == (mb, ma) else 'MISMATCH'}")
    out = dict(date=time.strftime("%Y-%m-%d %H:%M"), M_MAX=M_MAX, roots=res, levels=levels,
               targets=[f"J{d['J']}v{d['v']}" for d in targets], order_BOad=order_ad, order_BO=order_bo,
               bo_h2plus_check=chk, counts_M100=cnt)
    json.dump(out, open(os.path.join(HERE, "bo_map.json"), "w"), indent=1)
    P(f"# done in {time.time() - t0:.0f}s")
    open(os.path.join(HERE, "bo_map.txt"), "w").write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
