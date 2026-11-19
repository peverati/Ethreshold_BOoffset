"""Entanglement scan over the mass_scans grid.  *Last updated: 2026-11-13*

usage: python run.py MODE [M ...]
MODE: exact (base grid, largest basis) | fine (finer quadrature) | prevN (N-4 basis, base grid)
      | bo (BO product state, base grid) | dense (any M: ground state re-solved at M in the fixed basis
      (alphas, N) of the nearest scan mass DENSE_BASE, grid of that mass; via paper1/fidelity/fidcore)
Appends one row per M to results_<MODE>.csv and skips masses already present (resumable);
saves spectra to data/ent_<MODE>_M<M>.npz.
"""
import csv
import math
import os
import sys
import time

import numpy as np

import entcore as ec

MASSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.15,
          1.2, 1.5, 2.0, 3.0, 4.5, 6.0, 9.5256, 12.0, 15.0, 20.0, 29.639, 50.0, 100.0, 206.768283, 500.0,
          1000.0, 1836.15267343, 4000.0, 1.0e4]
COLS = ["M", "kappa", "S", "S_lin", "S_orient", "S_rad", "S_0", "w0", "w2", "w4", "w6", "w_lmax", "lam_max", "K",
        "E_N", "neg", "MI", "radial_dim", "radial_discard", "direct_norm", "sum_l", "time_s"]


def done(fn):
    if not os.path.exists(fn):
        return set()
    with open(fn) as fh:
        return {round(float(r["M"]), 6) for r in csv.DictReader(fh)}


DENSE_BASE = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5]
_PROB = {}


def dense_amp(M):
    sys.path.insert(0, os.path.join(ec.ROOT, "paper1", "fidelity"))
    import fidcore
    from threebody.solver import Result
    Mb = min(DENSE_BASE, key=lambda x: abs(math.log(x / M)))
    if Mb not in _PROB:
        _, b = ec.load_exact(Mb)
        pb = fidcore.Problem(b["N"], tuple(b["alphas"]), tuple(b["weights"]))
        _PROB.clear()
        _PROB[Mb] = (pb, pb.mats(b["N"]))
    pb, (bs, S, T01, T2, V) = _PROB[Mb]
    E, c, _ = fidcore.ground(S, T01, T2, V, M)
    res = Result(system=None, basis=bs, energies=np.array([E]), vectors=c[:, None])
    return ec.exact_amp(res), Mb


def one(M, mode):
    level = "fine" if mode == "fine" else "base"
    if mode == "dense":
        amp, Mb = dense_amp(M)
        grid = ec.grid_for(Mb, level)
    else:
        grid = ec.grid_for(M, level)
    if mode in ("exact", "fine"):
        res, _ = ec.load_exact(M)
        amp = ec.exact_amp(res)
    elif mode == "prevN":
        amp, _ = ec.exact_amp_prevN(M)
    elif mode == "bo":
        amp, _ = ec.bo_amp(M)
    A, direct = ec.partial_waves(amp, grid)
    sp = ec.spectrum(A, grid.ls)
    rm = ec.radial_measures(A, grid.ls)
    w = sp["w"]
    row = dict(M=M, kappa=M ** -0.25, S=sp["S"], S_lin=sp["S_lin"], S_orient=sp["S_orient"], S_rad=sp["S_rad"],
               S_0=sp["S_l"][0], w0=w[0], w2=w[1], w4=w[2], w6=w[3], w_lmax=w[-1], lam_max=sp["lam_max"], K=sp["K"],
               E_N=rm["E_N"], neg=rm["neg"], MI=rm["MI"], radial_dim=rm["k"] * rm["kp"],
               radial_discard=rm["discard"], direct_norm=direct.sum(), sum_l=sp["norm_sum"])
    np.savez(os.path.join(ec.DATA, f"ent_{mode}_M{M:g}.npz"), ls=grid.ls, w=w, S_l=sp["S_l"],
             lam_sorted=sp["lam_sorted"], lam_blocks=sp["lam_blocks"], grid=str(grid.desc))
    return row


if __name__ == "__main__":
    mode = sys.argv[1]
    masses = [float(x) for x in sys.argv[2:]] or MASSES
    fn = f"results_{mode}.csv"
    have = done(fn)
    new = not os.path.exists(fn)
    with open(fn, "a", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=COLS)
        if new:
            wr.writeheader()
        for M in masses:
            if round(M, 6) in have:
                continue
            t0 = time.time()
            try:
                row = one(M, mode)
            except Exception as e:  # noqa: BLE001
                print(f"M={M}: FAILED {e!r}", flush=True)
                continue
            row["time_s"] = time.time() - t0
            wr.writerow({k: (f"{v:.14g}" if isinstance(v, float) else v) for k, v in row.items()})
            fh.flush()
            print(f"{mode} M={M:g} S={row['S']:.8f} S_lin={row['S_lin']:.8f} Sor={row['S_orient']:.6f} "
                  f"Srad={row['S_rad']:.6f} E_N={row['E_N']:.6f} MI={row['MI']:.6f} norm={row['direct_norm']:.10f} "
                  f"t={row['time_s']:.0f}s", flush=True)
