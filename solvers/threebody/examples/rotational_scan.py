"""Mass scan of the lowest natural-parity L = 0, 1, 2 levels of {X+, X+, e-}.

For each M: L = 0 'S', L = 1 'A', L = 2 'S' (ground vibronic manifold, P = (-1)^L), and
L = 1 'S' (other exchange symmetry, for completeness). Scales optimized for root 0 at
N_opt, then a fixed-scale convergence sequence; two roots per symmetry.
Results are appended to <outdir>/scan_L.csv one (M, L, sym) row at a time (restartable:
rows already present are skipped).

Usage: python examples/rotational_scan.py [outdir] [M1 M2 ...]
"""
import csv
import json
import os
import sys
import time

import numpy as np

from threebody import symmetric_family, family_threshold
from threebody.rotational import convergence_table_L, optimize_scales

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                         "paper_data", "rotational")
MASSES = [float(x) for x in sys.argv[2:]] or [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0,
                                              206.768283, 500.0, 1000.0, 1836.15267343]
CASES = [(0, "S"), (1, "A"), (2, "S"), (1, "S")]


def settings(M):
    if M < 150:
        return dict(weights=(1, 1, 1), N_opt=12, Ns=[14, 18, 21])
    return dict(weights=(1, 1, 2), N_opt=16, Ns=[24, 28, 32])


FIELDS = ["M", "L", "sym", "E_th", "N_last", "size_last", "alpha_X", "alpha_e", "E0", "E0_ext", "E0_err", "E0_flag",
          "E1", "E1_ext", "E1_err", "E1_flag", "E0_seq", "E1_seq", "r01", "r01_2", "inv_r01", "virial", "time_s"]


def done_keys(path):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {(float(r["M"]), int(r["L"]), r["sym"]) for r in csv.DictReader(fh)}


def main():
    path = os.path.join(OUT, "scan_L.csv")
    have = done_keys(path)
    new = not os.path.exists(path)
    fh = open(path, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
    x0 = {}
    for M in MASSES:
        s = symmetric_family(M)
        st = settings(M)
        cur = {}
        for L, sym in CASES:
            if (M, L, sym) in have:
                continue
            if (L, sym) == (1, "S") and M > 10.5:
                continue  # 2p sigma_u-like, diffuse; not attempted at large M
            t0 = time.time()
            # scales: optimized for L=0 'S' and L=1 'A' (root 0); L=2 'S' reuses the L=1 'A' scales and
            # L=1 'S' the L=0 ones (fixed scales are still variational; convergence in N is checked)
            src = {(0, "S"): (0, "S"), (1, "A"): (1, "A"), (2, "S"): (1, "A"), (1, "S"): (0, "S")}[(L, sym)]
            if src == (L, sym) or src not in cur:
                guess = x0.get(src)
                al, e, n = optimize_scales(s, src[0], st["N_opt"], src[1], x0=guess, weights=st["weights"])
                x0[src] = al
                cur[src] = al
            al = cur[src]
            tab = convergence_table_L(s, L, sym, al, st["Ns"], nstates=2, weights=st["weights"], observables=True)
            last = tab["rows"][-1]
            ob = tab["observables"][0]
            row = dict(M=M, L=L, sym=sym, E_th=family_threshold(M), N_last=last["N"], size_last=last["size"],
                       alpha_X=al[0], alpha_e=al[2])
            for k in range(2):
                E, err, rho, flag = tab["extrap"][k]
                row.update({f"E{k}": repr(float(last["energies"][k])), f"E{k}_ext": repr(float(E)),
                            f"E{k}_err": f"{err:.2e}", f"E{k}_flag": flag,
                            f"E{k}_seq": json.dumps([float(r["energies"][k]) for r in tab["rows"]])})
            row.update(r01=ob["r01"], r01_2=ob["r01_2"], inv_r01=ob["inv_r01"], virial=ob["virial"],
                       time_s=round(time.time() - t0, 1))
            w.writerow(row)
            fh.flush()
            print(f"M={M:g} L={L}{sym}: E0={last['energies'][0]:.12f} (+/-{tab['extrap'][0][1]:.1e}) "
                  f"E1={last['energies'][1]:.10f} E_th={family_threshold(M):.10f} t={time.time() - t0:.0f}s",
                  flush=True)
    fh.close()


if __name__ == "__main__":
    main()
