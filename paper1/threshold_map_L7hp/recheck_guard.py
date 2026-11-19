"""Apply run 6's spurious-root guards to every finished crossing row of a threshold_map_L7hp target (2026-11-17).

For each row of data/<T>_crossings.csv with a finite M*, the basis is rebuilt (same code path as run_map_hp.py)
and checked:
  (c) crossing slope: gap(M*(1 - d)) + margin >= GUARD_BELOW and
      [gap(M*(1 + d)) - gap(M*(1 - d))] < -GUARD_SLOPE, d = 1e-3 (a genuine crossing; a spurious near-threshold
      level is flat or already bound below M*);
  conditioning: projected-S lam_min/lam_max < COND_FLAG (-1e-16) flags the basis (noise directions in S).
  (The ground-state check (b) of run 6 is for k >= 1 only: for k = 0 root 0 is the threshold state itself.)
Rows with a non-finite M* (no bracket) are checked at the bracket ends instead: a gap < 0 at the LOWER bracket end
(bound below every genuine crossing of a superset basis) marks the basis spurious.
Output: data/<T>_guard.csv (label, tol, margin, M_star, gap_lo, gap_hi, cond, flag, reason).
Usage (this folder): HP_NPROC=20 python recheck_guard.py J7v0"""
import csv
import json
import os
import sys

import numpy as np

import run_map_hp as R

D, GUARD_BELOW, GUARD_SLOPE, COND_FLAG = 1e-3, -1e-10, 1e-9, -1e-16


def parse(label):
    N = int(label.split("N=")[1].split()[0])
    nd = int(label.split("+ ")[1].split()[0]) if "diffuse" in label else 0
    return N, nd


def main(t):
    cfg = R.config(t)
    sc = json.load(open(os.path.join(R.DATA, f"{t}_scales.json")))
    mu = sc["M_opt"] / (sc["M_opt"] + 1)
    rows = list(csv.DictReader(open(os.path.join(R.DATA, f"{t}_crossings.csv"))))
    lo, hi = R.bracket(cfg)
    out = []
    for lab in dict.fromkeys(r["label"] for r in rows):
        N, nd = parse(lab)
        b = R.make_basis(cfg, sc["scales"], N, nd, mu)
        for r in [r for r in rows if r["label"] == lab]:
            tol, m, Ms = float(r["tol"]), float(r["margin"]), float(r["M_star"])
            b._proj = None
            pr = b.project(tol)
            cond = pr["lam_min"] / pr["lam_max"]
            why = []
            if np.isfinite(Ms):
                glo = R.gaps(cfg, b, Ms * (1 - D), tol)[cfg["k"]] + m
                ghi = R.gaps(cfg, b, Ms * (1 + D), tol)[cfg["k"]] + m
                if glo < GUARD_BELOW or ghi - glo > -GUARD_SLOPE:
                    why.append(f"crossing slope: gap(M*(1-{D:g}))={glo:+.3e}, gap(M*(1+{D:g}))={ghi:+.3e}")
            else:
                glo = R.gaps(cfg, b, lo, tol)[cfg["k"]] + m
                ghi = R.gaps(cfg, b, hi, tol)[cfg["k"]] + m
                if glo < 0:
                    why.append(f"bound at the lower bracket end M={lo:.4f} (gap {glo:+.3e}): spurious level")
            if cond < COND_FLAG:
                why.append(f"conditioning: lam_min/lam_max={cond:+.1e} < {COND_FLAG:g}")
            out.append(dict(label=lab, tol=r["tol"], margin=r["margin"], M_star=r["M_star"], gap_lo=repr(glo),
                            gap_hi=repr(ghi), cond=repr(cond), flag=int(bool(why)), reason="; ".join(why)))
            R.log(t, f"guard recheck {lab} tol={tol:g} margin={m:g} M*={Ms:.10f}: "
                     f"{'FLAGGED: ' + '; '.join(why) if why else 'ok'} (gap_lo {glo:+.3e}, gap_hi {ghi:+.3e}, cond {cond:+.1e})")
        del b
    with open(os.path.join(R.DATA, f"{t}_guard.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main(sys.argv[1])
