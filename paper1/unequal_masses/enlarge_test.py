"""Diffuse-basis enlargement test for the L = 0 unequal-mass targets (written 2026-11-16, Mac Studio).

Why: the r = 1 check R1J0v1 gave 9.5873 (production basis) vs the symmetric exact 9.525634, and the production
"diffuse enlarged" variant (n_c + 2, n_d + 20 on all diffuse blocks) dropped it to 9.5590.  Hypotheses:
(H1) at r = 1 the "closed" channel (B far from A e) is degenerate with the open one and needs the same diffuse
ladder ("mirror"); (H2) the atom polarization (degree n_c in the atom coordinates of the diffuse blocks) is too
low ("pol").  Variants (compact N = 16 with the stage-1 scales and the standard open blocks always present):
  M    mirror: closed channel also gets (4,16,0.15), (4,30,0.03), (4,50,alpha_s=1e-4)
  P6   n_c 4 -> 6 on every diffuse block
  P8   n_c 4 -> 8 on every diffuse block
  MP6  mirror + n_c 6
  MP8  mirror + n_c 8
  MB   mirror + (n_c + 2, n_d + 20) on every diffuse block
  MD20 / MD40  mirror + n_d + 20 / + 40 on every diffuse block
  MN14 / MN18  mirror with compact N = 14 / 18 (compact-size sequence of the mirror basis)
Name suffix "@m" (e.g. MP6@m): crossing at margin 1e-11 instead of 0 (for the largest bases the margin-0 gap
below M_c is at the 1e-13 noise level, so a margin-0 search can stop on noise; MP6 at r = 1 did).  Env ENL_LO
overrides the lower bracket.  The shift column is relative to the production margin-0 value in all cases.
Each variant: one crossing (tol 1e-12) with the run_exact.py guard (c); results in
enlarge_<target>.csv and enlarge_<target>_log.txt (resumable: done variants are skipped).
Run from the project root:
    THREEBODY_BLAS_THREADS=5 python paper1/unequal_masses/enlarge_test.py R1J0v1 M P6 MP6 P8 MB
"""
import csv
import gc
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402

import run_exact as R  # noqa: E402

ALPHA = 1e-4
MIRROR = ((4, 16, 0.15), (4, 30, 0.03), (4, 50, ALPHA))        # extra closed-channel sets (mirror of open)
VARIANTS = {
    "M": dict(mirror=True, dc=0, dd=0),
    "P6": dict(mirror=False, dc=2, dd=0),
    "P8": dict(mirror=False, dc=4, dd=0),
    "MP6": dict(mirror=True, dc=2, dd=0),
    "MP8": dict(mirror=True, dc=4, dd=0),
    "MB": dict(mirror=True, dc=2, dd=20),
    "MD20": dict(mirror=True, dc=0, dd=20),
    "MD40": dict(mirror=True, dc=0, dd=40),
    "MN14": dict(mirror=True, dc=0, dd=0, N=14),
    "MN18": dict(mirror=True, dc=0, dd=0, N=18),
}
MARGIN = 1e-11   # variant name suffix "@m": crossing of E_k = E_th - MARGIN (noise-robust; the certified bound)
COLS = ["target", "variant", "size", "kept", "M_star", "M_std", "shift", "spurious", "reason", "seconds", "finished"]


def blocks(t, sc, mirror, dc, dd, N=16):
    cfg = R.TARGETS[t]
    r, Mref = cfg["r"], sc["M_opt"]
    muA, muB = Mref / (Mref + 1), r * Mref / (r * Mref + 1)
    bl = [R.block(N, tuple(a), name=f"N{N}s{i}") for i, a in enumerate(sc["scales"])]
    bl.append(R.open_block(ALPHA, muB, 4 + dc, cfg["nd"] + dd))
    bl += [R.open_block(a, muB, nc + dc, nd + dd) for nc, nd, a in R.OPEN_FIXED]
    closed = list(R.CLOSED) + (list(MIRROR) if mirror else [])
    bl += [R.closed_block(a, muA, nc + dc, nd + dd) for nc, nd, a in closed]
    return bl


def std_value(t):
    rows = list(csv.DictReader(open(os.path.join(HERE, f"{t}_crossings.csv"))))
    lab = f"compact N=16 + diffuse alpha_s={ALPHA:g}"
    return [float(r["M_star"]) for r in rows if r["label"] == lab and r["tol"] == repr(1e-12)
            and r["margin"] == repr(0.0)][0]


def main():
    t, names = sys.argv[1], sys.argv[2:]
    R._logs[t] = open(os.path.join(HERE, f"enlarge_{t}_log.txt"), "a")     # keep production logs untouched
    path = os.path.join(HERE, f"enlarge_{t}.csv")
    done = {r["variant"] for r in csv.DictReader(open(path))} if os.path.exists(path) else set()
    sc = R.stage_scales(t)
    lo, _, _ = R.brackets(t)
    lo = float(os.environ.get("ENL_LO", lo))
    Mstd = std_value(t)
    for name in names:
        if name in done:
            continue
        v = VARIANTS[name.split("@")[0]]
        m = MARGIN if name.endswith("@m") else 0.0
        t0 = time.perf_counter()
        bl = blocks(t, sc, **v)
        R.log(t, f"enlarge {name} ({v}): size {sum(b.size for b in bl)}")
        b = R.LeanBasis(0, bl, R.TARGETS[t]["r"])
        try:
            Ms = R.crossing(t, b, lo, Mstd * (1 + 1e-4), margin=m)
        except ValueError as exc:
            R.log(t, f"  enlarge {name}: {exc}")
            Ms = float("nan")
        why = R.guard(t, b, Ms, 1e-12, v.get("N", 16), margin=m)
        row = dict(target=t, variant=name, size=b.size, kept=b.proj[1e-12]["kept"],
                   M_star=float("nan") if why else Ms, M_std=Mstd, shift=(float("nan") if why else Ms - Mstd),
                   spurious=int(bool(why)), reason=(f"M*={Ms:.10f}; {why}" if why else ""),
                   seconds=round(time.perf_counter() - t0), finished=time.strftime("%Y-%m-%d %H:%M:%S"))
        new = not os.path.exists(path)
        with open(path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS)
            if new:
                w.writeheader()
            w.writerow({k: (repr(x) if isinstance(x, float) else x) for k, x in row.items()})
        R.log(t, f"enlarge {name}: size={b.size} kept={row['kept']} M*={Ms:.10f} (std {Mstd:.10f}, shift "
                 f"{Ms - Mstd:+.3e}){' REJECTED: ' + why if why else ''} ({row['seconds']}s)")
        del b
        gc.collect()
    R.log(t, "enlarge DONE")


if __name__ == "__main__":
    np.seterr(all="ignore")
    main()
