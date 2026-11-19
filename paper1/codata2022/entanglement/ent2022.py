"""Electron-nuclear entanglement at the CODATA-2022 physical masses (H2+, Mu2+).

*Last updated: 2026-11-17 17:48*

Uses paper1/entanglement (entcore, run.one, bodyframe.body_entropy) unchanged, but reads the exact
wave functions from ../data (the CODATA-2022 scan) and writes spectra to ./data and the CSVs here,
so that the CODATA-2018 results in paper1/entanglement are not touched.
Usage (from anywhere):  python ent2022.py      (modes exact, prevN, bo; body frame for H2+)
"""
import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "paper1", "entanglement"))
sys.path.insert(0, os.path.join(ROOT, "paper1", "codata2022"))

import entcore as ec  # noqa: E402
from common import M_PROTON, M_MUON  # noqa: E402

ec.SCANDATA = os.path.join(ROOT, "paper1", "codata2022", "data")
ec.DATA = os.path.join(HERE, "data")
os.makedirs(ec.DATA, exist_ok=True)

import run  # noqa: E402

for mode in ("exact", "prevN", "bo"):
    fn = os.path.join(HERE, f"results_{mode}.csv")
    have = run.done(fn)
    new = not os.path.exists(fn)
    with open(fn, "a", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=run.COLS)
        if new:
            wr.writeheader()
        for M in (M_MUON, M_PROTON):
            if round(M, 6) in have:
                continue
            t0 = time.time()
            row = run.one(M, mode)
            row["time_s"] = time.time() - t0
            wr.writerow({k: (f"{v:.14g}" if isinstance(v, float) else v) for k, v in row.items()})
            fh.flush()
            print(f"{mode} M={M:.12g} S={row['S']:.10f} Sor={row['S_orient']:.8f} Srad={row['S_rad']:.8f} "
                  f"t={row['time_s']:.0f}s", flush=True)

# body-frame entropy (Sanz-Vicario comparison) at the CODATA-2022 proton mass
sys.argv = [sys.argv[0], repr(M_PROTON)]
ec.HERE = HERE  # bodyframe.py writes results_bodyframe.csv into ec.HERE
import bodyframe  # noqa: E402,F401  (runs on import, as in paper1/entanglement)
