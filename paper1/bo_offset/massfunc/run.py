"""Resumable production driver for m_v(R), m_r(R) (see README.md, THEORY.md).

Run from the project root:
    python3 paper1/bo_offset/massfunc/run.py [--fine] [--dry] [--grid FILE] [--R 1 2 3]
Default: ~650 R points in [0.02, 60] bohr, two basis levels per point (0 and 1; value = level 1,
convergence = |level1 - level0|) -> results.csv.  --fine: ~2600 points, levels 1 and 2 -> results_fine.csv.
One CSV row is appended (and flushed) per R; R values already in the output file are skipped.
"""
import argparse
import csv
import os
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
os.environ.setdefault("MASSFUNC_THREADS", "4")
import numpy as np

COLS = ["R", "m_v", "m_r", "dm_v", "dm_r", "m_v_lo", "m_r_lo", "A_par", "A_perp", "E0_exact", "dE0_basis",
        "fpar_ratio", "fperp_ratio", "E1_sigma_g", "E_pi_g", "n_xi", "n_eta", "n_xi_lo", "n_eta_lo",
        "seconds", "finished"]


def grid(fine=False):
    s = 0.25 if fine else 1.0
    segs = [(0.02, 0.5, 0.02 * s), (0.5, 10.0, 0.02 * s), (10.0, 20.0, 0.1 * s), (20.0, 40.0, 0.5 * s),
            (40.0, 60.0 + 1e-9, 2.0 * s)]
    pts = []
    for a, b, h in segs:
        pts.extend(np.arange(a, b - 1e-9, h))
    pts.append(60.0)
    return np.unique(np.round(pts, 6))


def done_R(path):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {round(float(r["R"]), 6) for r in csv.DictReader(fh)}


def estimate(Rs, levels):
    """Wall-time estimate from a one-point calibration (cost ~ basis size^3)."""
    import massfunc as mf
    t0 = time.time(); mf.mass_functions(10.0, level=0); t = time.time() - t0
    n10 = np.prod(mf.default_sizes(10.0, 0))
    tot = 0.0
    for R in Rs:
        for lv in levels:
            n = np.prod(mf.default_sizes(R, lv))
            tot += t * max((n / n10) ** 3, 0.3)
    return tot


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine", action="store_true", help="denser grid, larger bases (levels 1,2)")
    ap.add_argument("--dry", action="store_true", help="print plan and estimated time only")
    ap.add_argument("--grid", help="file with one R per line (overrides built-in grid)")
    ap.add_argument("--R", type=float, nargs="*", help="explicit R values (overrides grid)")
    ap.add_argument("--out", help="output CSV (default results.csv / results_fine.csv in this folder)")
    a = ap.parse_args(argv)
    levels = (1, 2) if a.fine else (0, 1)
    out = a.out or os.path.join(here, "results_fine.csv" if a.fine else "results.csv")
    if a.R:
        Rs = np.array(a.R)
    elif a.grid:
        Rs = np.loadtxt(a.grid, ndmin=1)
    else:
        Rs = grid(a.fine)
    have = done_R(out)
    todo = [R for R in Rs if round(float(R), 6) not in have]
    print(f"output {out}\n{len(Rs)} R points in [{Rs.min()}, {Rs.max()}], {len(Rs) - len(todo)} already done, "
          f"{len(todo)} to do; basis levels {levels}; threads {os.environ['MASSFUNC_THREADS']}", flush=True)
    if a.dry:
        est = estimate(todo, levels)
        print(f"estimated wall time {est / 60:.1f} min (one-point calibration on this machine; load-dependent)")
        return
    import massfunc as mf
    new = not os.path.exists(out)
    fh = open(out, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=COLS)
    if new:
        w.writeheader(); fh.flush()
    t_start = time.time()
    for k, R in enumerate(todo):
        t0 = time.time()
        lo = mf.mass_functions(float(R), level=levels[0])
        hi = mf.mass_functions(float(R), level=levels[1])
        row = dict(R=f"{R:.6f}", m_v=f"{hi['m_v']:.14f}", m_r=f"{hi['m_r']:.14f}",
                   dm_v=f"{abs(hi['m_v'] - lo['m_v']):.3e}", dm_r=f"{abs(hi['m_r'] - lo['m_r']):.3e}",
                   m_v_lo=f"{lo['m_v']:.14f}", m_r_lo=f"{lo['m_r']:.14f}",
                   A_par=f"{hi['A_par']:.14f}", A_perp=f"{hi['A_perp']:.14f}",
                   E0_exact=f"{hi['E0_exact']:.14f}", dE0_basis=f"{hi['E0_basis'] - hi['E0_exact']:.3e}",
                   fpar_ratio="", fperp_ratio="",
                   E1_sigma_g=f"{hi['E1_sigma_g']:.12f}", E_pi_g=f"{hi['E_pi_g']:.12f}",
                   n_xi=hi["n_xi"], n_eta=hi["n_eta"], n_xi_lo=lo["n_xi"], n_eta_lo=lo["n_eta"],
                   seconds=f"{time.time() - t0:.2f}", finished=time.strftime("%Y-%m-%d %H:%M:%S"))
        if k % 25 == 0 or R in (2.0, 30.0, 60.0):      # occasional completeness check against D(R)
            dr, dp = mf.reference_norms(float(R))
            row["fpar_ratio"] = f"{hi['fpar_norm_in_basis'] / dr:.8f}"
            row["fperp_ratio"] = f"{hi['fperp_norm_in_basis'] / dp:.8f}"
        w.writerow(row); fh.flush(); os.fsync(fh.fileno())
        if k % 20 == 0 or k == len(todo) - 1:
            el = time.time() - t_start
            print(f"[{k + 1}/{len(todo)}] R={R:.4f} m_v={hi['m_v']:.10f} m_r={hi['m_r']:.10f} "
                  f"dm={float(row['dm_v']):.1e}/{float(row['dm_r']):.1e} elapsed {el / 60:.1f} min, "
                  f"ETA {el / (k + 1) * (len(todo) - k - 1) / 60:.1f} min", flush=True)
    fh.close()
    print("done", out, flush=True)


if __name__ == "__main__":
    main()
