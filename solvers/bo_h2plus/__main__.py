"""Command line interface:  python -m bo_h2plus <command> ...   (run from the solvers/ folder)

  electronic --R 2.0 [--state 1sg|2pu] [--mu 1.0]
  equilibrium [--M 1836.15267343]
  levels --M 1836.15267343 [--J 0] [--adiabatic] [--out levels.csv] [--density dens.csv]
  scan --M 1 10 100 1836.15267343 1e4 --out scan.csv [--adiabatic]
  scan --logspace 0 5 21 --out scan.csv
  critical [--nmax 8]
"""
import argparse
import csv
import sys

import numpy as np


def main(argv=None):
    from . import electronic as el, potential as P
    from .nuclear import solve_levels, count_levels, thresholds
    from .scan import mass_scan, critical_masses, HARTREE_CM, m_eff

    ap = argparse.ArgumentParser(prog="bo_h2plus", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("electronic", help="clamped-nuclei energy at R")
    a.add_argument("--R", type=float, nargs="+", required=True)
    a.add_argument("--state", default="1sg", choices=["1sg", "2pu"])
    a.add_argument("--mu", type=float, default=1.0, help="electron (reduced) mass")
    a = sub.add_parser("equilibrium", help="R_e, E_min, k (for mu_e of mass M)")
    a.add_argument("--M", type=float, default=None)
    a = sub.add_parser("levels", help="all bound levels for M, J")
    a.add_argument("--M", type=float, required=True)
    a.add_argument("--J", type=int, default=0)
    a.add_argument("--adiabatic", action="store_true")
    a.add_argument("--out")
    a.add_argument("--density", help="CSV of P_v(R)=|u_v(R)|^2 on the DVR grid")
    a = sub.add_parser("scan", help="mass scan to CSV")
    a.add_argument("--M", type=float, nargs="*")
    a.add_argument("--logspace", type=float, nargs=3, metavar=("LOG10_MIN", "LOG10_MAX", "N"))
    a.add_argument("--out", required=True)
    a.add_argument("--adiabatic", action="store_true")
    a = sub.add_parser("critical", help="BO critical masses")
    a.add_argument("--nmax", type=int, default=8)
    a.add_argument("--out")
    args = ap.parse_args(argv)

    if args.cmd == "electronic":
        for R in args.R:
            Ee = el.electronic_energy_mu(R, args.mu, args.state)
            print(f"R={R:.10g} mu={args.mu:.12g} state={args.state}  E_el={Ee:.15f}  "
                  f"E_el+1/R={Ee + 1 / R:.15f}")
    elif args.cmd == "equilibrium":
        mu = 1.0 if args.M is None else P.mu_e(args.M)
        eq = el.equilibrium(mu=mu)
        print(f"mu_e={mu:.15g}  R_e={eq.Re:.12f}  E_min={eq.Emin:.15f}  k={eq.k:.12f}")
    elif args.cmd == "levels":
        s = solve_levels(args.M, args.J, adiabatic=args.adiabatic)
        thr_bo, thr_ex = thresholds(args.M, args.adiabatic)
        print(f"# M={args.M} J={args.J} adiabatic={args.adiabatic}  BO threshold={thr_bo:.15f} "
              f"exact 2-body threshold={thr_ex:.15f}")
        print(f"# bound levels: DVR {s.n_dvr}, node count {count_levels(args.M, args.J, args.adiabatic)}"
              f"  (grid points {s.grid_info['npts']}, R_max {s.grid_info['Rmax']:.4g})")
        cols = ["v", "J", "E", "binding", "err_E", "R_mean", "dR", "lindemann", "inv_R", "inv_R2",
                "R2_mean", "T_kin", "T_rad"]
        print(" ".join(f"{c:>18s}" for c in cols))
        for l in s.levels:
            print(" ".join(f"{getattr(l, c):18.12g}" if c not in ("v", "J") else f"{getattr(l, c):18d}"
                           for c in cols))
        if args.out:
            with open(args.out, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(cols)
                for l in s.levels:
                    w.writerow([repr(getattr(l, c)) if isinstance(getattr(l, c), float) else getattr(l, c)
                                for c in cols])
        if args.density:
            np.savetxt(args.density, np.column_stack([s.R, s.P.T]), delimiter=",",
                       header="R," + ",".join(f"P_v{v}" for v in range(len(s.levels))), comments="")
    elif args.cmd == "scan":
        masses = list(args.M or [])
        if args.logspace:
            lo, hi, n = args.logspace
            masses += list(np.logspace(lo, hi, int(n)))
        mass_scan(masses, out=args.out, adiabatic=args.adiabatic, verbose=True)
        print("wrote", args.out)
    elif args.cmd == "critical":
        res = critical_masses(args.nmax)
        print("# BO prediction (no adiabatic correction): mass M at which level v appears at -mu_e/2")
        print(f"{'n':>3} {'v':>3} {'M_c':>22} {'m_eff=(2M+1)/4':>20} {'est_abs_err':>12}")
        for d in res:
            print(f"{d['n']:3d} {d['v_new']:3d} {d['M_c']:22.12f} {d['m_eff']:20.12f} {d['est_abs_err']:12.2e}")
        if args.out:
            with open(args.out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(res[0].keys()))
                w.writeheader()
                w.writerows(res)
    return 0


if __name__ == "__main__":
    sys.exit(main())
