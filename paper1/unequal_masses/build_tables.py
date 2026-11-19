"""Build the (beta-independent) tables used by the theory: g/u coupling functions and mass matrices.

Written 2026-11-16.  Run from the project root (about 5-10 min, < 1 GB):
    python3 paper1/unequal_masses/build_tables.py
Outputs: solvers/bo_h2plus/_cache/gu_v1.npz (g/u coupling functions, via bo_h2plus.unequal.gu_table) and
paper1/unequal_masses/massfunc_gu.csv (mass-matrix coefficients; level 0, with level-1 values every 5th
point for the basis-convergence column d_max).  Resumable: R values already in the CSV are skipped.
"""
import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402

import massfunc_gu as MG  # noqa: E402
from bo_h2plus import unequal as UQ  # noqa: E402

OUT = os.path.join(HERE, "massfunc_gu.csv")


def main():
    t0 = time.time()
    tab = UQ.gu_table(verbose=True)
    print(f"gu table: {len(tab['R'])} points ({time.time() - t0:.0f}s)", flush=True)
    done = set()
    if os.path.exists(OUT):
        done = {round(float(r["R"]), 8) for r in csv.DictReader(open(OUT))}
    cols = ["R"] + MG.KEYS + ["d_max", "E_Sg1", "E_Su1", "E_Pg0", "E_Pu0", "e_g", "e_u", "seconds"]
    new = not os.path.exists(OUT)
    with open(OUT, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        if new:
            w.writeheader()
        for i, R in enumerate(MG.default_grid()):
            if round(float(R), 8) in done:
                continue
            q = MG.mass_matrices(R, level=0)
            q["d_max"] = float("nan")
            if i % 5 == 0:
                q1 = MG.mass_matrices(R, level=1)
                q["d_max"] = max(abs(q1[k] - q[k]) for k in MG.KEYS)
            w.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in q.items()})
            fh.flush()
            if i % 10 == 0:
                print(f"  R={R:.3f} mv_gg0={q['mv_gg0']:.6f} mv_gu1={q['mv_gu1']:.6f} d_max={q['d_max']:.1e}", flush=True)
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
