"""Command line interface.

Examples
--------
python -m threebody solve --M 1 --N 20 --nstates 4
python -m threebody solve --masses 1 1 inf --charges -1 -1 2 --identical --N 20
python -m threebody converge --M 1836.15267343 --N 12 16 20 24 --nstates 3
python -m threebody scan --M 1 2 5 10 --N 12 16 20 --nstates 4 --out scan.csv --save-dir scan_out
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

from . import (ThreeBodySystem, PerimetricBasis, solve, symmetric_family, optimize_alphas,
               convergence_table)
from .convergence import format_table, default_symmetry
from .scan import scan_masses, write_csv, write_gr


def _system(a):
    if a.masses is not None:
        masses = [float(x) for x in a.masses]
        charges = [float(x) for x in a.charges]
        return ThreeBodySystem(masses, charges, identical=(0, 1) if a.identical else None)
    return symmetric_family(a.M[0])


def _sym(a, system):
    if a.symmetry == "auto":
        return default_symmetry(system)
    return None if a.symmetry == "none" else a.symmetry


def main(argv=None):
    p = argparse.ArgumentParser(prog="threebody", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("solve", "converge", "scan"):
        q = sub.add_parser(name)
        q.add_argument("--M", type=float, nargs="+", default=[1.0], help="X mass(es) for the {X+,X+,e-} family")
        q.add_argument("--masses", nargs=3, default=None, help="general system: three masses (inf allowed)")
        q.add_argument("--charges", nargs=3, default=None)
        q.add_argument("--identical", action="store_true", help="particles 0 and 1 identical")
        q.add_argument("--symmetry", default="auto", choices=["auto", "S", "A", "none"])
        q.add_argument("--N", type=int, nargs="+", default=[16])
        q.add_argument("--N-opt", type=int, default=12, help="basis size used to optimize the scales")
        q.add_argument("--alphas", type=float, nargs=3, default=None, help="fixed scales (skip optimization)")
        q.add_argument("--nstates", type=int, default=4)
        q.add_argument("--opt-state", type=int, default=0, help="root whose energy is minimized to choose the scales")
        q.add_argument("--weights", type=int, nargs=3, default=[1, 1, 1],
                       help="truncation w0 a0 + w1 a1 + w2 a2 <= N (default total degree)")
        q.add_argument("--grid", default=None, help="file with r values (one per line) for g(r)")
        q.add_argument("--out", default=None)
        q.add_argument("--save-dir", default=None)
    a = p.parse_args(argv)

    grid = np.loadtxt(a.grid, ndmin=1) if a.grid else None

    if a.cmd == "scan":
        rows = scan_masses(a.M, a.N, a.nstates, a.N_opt, symmetry="S" if a.symmetry == "auto" else _sym(a, None),
                           save_dir=a.save_dir, alphas=a.alphas, grid=grid, opt_state=a.opt_state,
                           weights=tuple(a.weights))
        write_csv(a.out or "scan.csv", rows)
        print(f"wrote {a.out or 'scan.csv'}")
        return 0

    system = _system(a)
    sym = _sym(a, system)
    if a.alphas is not None:
        al = tuple(a.alphas)
    else:
        al, e, nev = optimize_alphas(system, a.N_opt, sym, state=a.opt_state, weights=tuple(a.weights))
        print(f"optimized alphas at N={a.N_opt}: {np.round(al, 8).tolist()}  ({nev} evaluations)")
    if a.cmd == "converge":
        tab = convergence_table(system, al, a.N, sym, a.nstates, tuple(a.weights))
        print(format_table(tab))
        res = tab["last"]
    else:
        res = solve(system, PerimetricBasis(a.N[-1], al, sym, tuple(a.weights)), a.nstates)
        print(f"basis size {res.basis.size}; threshold {res.threshold:.12f}")
        for k, e in enumerate(res.energies):
            print(f"  E{k} = {e:.14f}   E - E_th = {e - res.threshold:+.3e}")
    print(f"variational roots below threshold (=> at least this many bound states): {res.n_below_threshold()}")
    obs = res.observables(0)
    for k, v in obs.items():
        print(f"  {k:12s} {v:.12g}")
    if a.save_dir:
        import os
        os.makedirs(a.save_dir, exist_ok=True)
        res.save(os.path.join(a.save_dir, "result.npz"))
        if grid is not None:
            write_gr(os.path.join(a.save_dir, "gr.txt"), res, grid)
    return 0


if __name__ == "__main__":
    sys.exit(main())
