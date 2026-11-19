"""Combine the exact L = 0, 1, 2 scan (scan_L.csv) with the BO / BO+adiabatic levels of bo_h2plus.

Writes <outdir>/rotational_table.csv and rotational_table.txt.
Run from solvers/threebody with solvers/ on the path (done below).
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))  # solvers/
from bo_h2plus import solve_levels  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "..", "paper1", "rotational")
CM = 219474.6313632


def bo_levels(M, adiabatic):
    out = {}
    for J in (0, 1, 2):
        try:
            ls = solve_levels(M, J, adiabatic=adiabatic)
            out[J] = ([lv.E for lv in ls.levels], ls.threshold)
        except Exception as exc:  # noqa: BLE001
            out[J] = ([], np.nan)
    return out


def main():
    rows = list(csv.DictReader(open(os.path.join(OUT, "scan_L.csv"))))
    by = {}
    for r in rows:
        by.setdefault(float(r["M"]), {})[(int(r["L"]), r["sym"])] = r
    table = []
    for M in sorted(by):
        d = by[M]
        if not all(k in d for k in ((0, "S"), (1, "A"), (2, "S"))):
            continue
        g = lambda k, i=0: float(d[k][f"E{i}"])
        err = lambda k, i=0: float(d[k][f"E{i}_err"])
        Eth = float(d[(0, "S")]["E_th"])
        E00, E10 = g((0, "S")), g((0, "S"), 1)
        E01, E02 = g((1, "A")), g((2, "S"))
        E11, E12 = g((1, "A"), 1), g((2, "S"), 1)
        row = dict(M=M, E_th=Eth, E00=E00, E01=E01, E02=E02, E10=E10, E11=E11, E12=E12,
                   err00=err((0, "S")), err01=err((1, "A")), err02=err((2, "S")),
                   bound01=E01 < Eth - 1e-10, bound02=E02 < Eth - 1e-10, bound10=E10 < Eth - 1e-10,
                   E01_S=g((1, "S")) if (1, "S") in d else np.nan)
        row["d1"] = E01 - E00
        row["d2"] = E02 - E00
        row["ratio"] = row["d2"] / row["d1"]
        row["dvib"] = E10 - E00
        row["rigid_B_from_inv_r2"] = np.nan
        for tag, ad in (("BO", False), ("BOad", True)):
            lv = bo_levels(M, ad)
            e0 = lv[0][0]
            def lev(J, v):
                return lv[J][0][v] if len(lv[J][0]) > v else np.nan
            row[f"{tag}_th"] = lv[0][1]
            row[f"{tag}_E00"] = lev(0, 0)
            row[f"{tag}_d1"] = lev(1, 0) - lev(0, 0)
            row[f"{tag}_d2"] = lev(2, 0) - lev(0, 0)
            row[f"{tag}_ratio"] = row[f"{tag}_d2"] / row[f"{tag}_d1"]
            row[f"{tag}_dvib"] = lev(0, 1) - lev(0, 0)
        table.append(row)
    keys = list(table[0].keys())
    with open(os.path.join(OUT, "rotational_table.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in table:
            w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    with open(os.path.join(OUT, "rotational_table.txt"), "w") as fh:
        fh.write("Exact (non-BO) lowest natural-parity levels of {X+,X+,e-}; L=0,2 'S', L=1 'A' (P=(-1)^L).\n")
        fh.write("d1 = E(L=1)-E(L=0), d2 = E(L=2)-E(L=0) [hartree]; ratio = d2/d1 (rigid rotor 3); dvib = E(v=1,L=0)-E(v=0,L=0).\n")
        fh.write("bound: lowest root below E_th = -M/(2(M+1)) (variational proof); '-' = not bound in this basis.\n")
        fh.write("M*d1: rigid rotor d1 = 2B, B = <1/R^2>/(2 mu_N) = <1/R^2>/M, so M*d1 -> 2<1/R^2> (constant) at large M.\n")
        fh.write("Levels marked '-' are NOT bound (root above E_th): their 'spacings' are pseudo-continuum numbers, not physical.\n")
        fh.write("bv: v=1, L=0 root below E_th in this (ground-state-optimized) basis; at M=10 it is bound but needs root-1 scales (README sec. 6). BO/BOad: bo_h2plus Born-Oppenheimer (+adiabatic) J-spacings of v=0 (nan: J level not bound in BO).\n\n")
        hdr = (f"{'M':>10} {'E(L=0)':>16} {'E(L=1)':>16} {'E(L=2)':>16} {'b1':>2} {'b2':>2} {'bv':>2} {'d1':>12} {'d2':>12} "
               f"{'ratio':>7} {'M*d1':>9} {'dvib':>11} {'BO d1':>12} {'BOad d1':>12} {'BO ratio':>8} {'d1/BOad':>8}")
        fh.write(hdr + "\n")
        for r in table:
            fh.write(f"{r['M']:10.4f} {r['E00']:16.12f} {r['E01']:16.12f} {r['E02']:16.12f} "
                     f"{'y' if r['bound01'] else '-':>2} {'y' if r['bound02'] else '-':>2} {'y' if r['bound10'] else '-':>2} {r['d1']:12.5e} {r['d2']:12.5e} "
                     f"{r['ratio']:7.4f} {r['M'] * r['d1']:9.5f} {r['dvib']:11.4e} {r['BO_d1']:12.5e} {r['BOad_d1']:12.5e} "
                     f"{r['BO_ratio']:8.4f} {r['d1'] / r['BOad_d1']:8.5f}\n")
    print(open(os.path.join(OUT, "rotational_table.txt")).read())


if __name__ == "__main__":
    main()
