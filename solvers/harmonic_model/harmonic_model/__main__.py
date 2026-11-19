"""CSV of all closed-form quantities over a list of masses.

    python -m harmonic_model --model A --masses 0.1 1 10 100 1836 [--k 1 --kp 0] [-o out.csv]
    python -m harmonic_model --model B --g 1 --masses 0.1 1 10 100 1836
    python -m harmonic_model --model twobody --masses 0.1 1 1836
"""
import argparse
import csv
import sys

import numpy as np

from . import model_a as A, model_b as B, twobody as T


def rows(model, masses, k, kp, g):
    for M in masses:
        if model == "A":
            d = A.summary(M, k, kp)
            yield {key: float(v) for key, v in d.items()}
        elif model == "B":
            yield B.summary(M, g, k, kp)
        else:
            yield dict(M=M, E_harm=float(T.harmonic_energy(M, k)), chi_harm=float(T.harmonic_chi_fixed(M)),
                       E_coul=float(T.coulomb_energy(M)), chi_coul=float(T.coulomb_chi_fixed(M)),
                       chi_inv=0.0)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="harmonic_model", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", choices=["A", "B", "twobody"], default="B")
    ap.add_argument("--masses", type=float, nargs="+", default=[0.1, 1, 10, 100, 1836])
    ap.add_argument("--k", type=float, default=1.0)
    ap.add_argument("--kp", type=float, default=0.0)
    ap.add_argument("--g", type=float, default=1.0)
    ap.add_argument("-o", "--out", default="-")
    a = ap.parse_args(argv)
    data = list(rows(a.model, a.masses, a.k, a.kp, a.g))
    fh = sys.stdout if a.out == "-" else open(a.out, "w", newline="")
    w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
    w.writeheader()
    for r in data:
        w.writerow({kk: (f"{v:.12g}" if isinstance(v, float) else v) for kk, v in r.items()})
    if fh is not sys.stdout:
        fh.close()
    if a.model == "B":
        print(f"# M_c (CM-density curvature sign change) for g={a.g}: "
              f"{B.critical_mass_cm(a.g, a.k, a.kp):.10g}", file=sys.stderr)


if __name__ == "__main__":
    main()
