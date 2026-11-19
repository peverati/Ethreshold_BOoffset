"""Is run 5 (double precision, tol 1e-12, single-pass canonical orthogonalization) affected by the orthonormality
problem found for L = 7 at tol <= 1e-14?  Written 2026-11-17.

For representative run-5 targets the standard largest basis (compact N_D + all diffuse sets, stage-1 scales from
paper1/threshold_map/data) is rebuilt with run 5's own code (double-precision CoupledBasis; term lists from this
folder's private cache copy, identical to run 5's, so run 5's cache is never touched).  Reported: single-pass
max|X^T S X - I| at tol 1e-12, the second-pass eigenvalue range, and the crossing M* with one and with two passes
(same bracket and Brent settings as run_map.py).  Usage (this folder): python check_run5_orthog.py J2v1 J5v0 J6v1"""
import csv
import json
import os
import sys

import numpy as np

import tb_hp_local
sys.modules["tb_local"] = tb_hp_local          # run_map imports tb_local: redirect to the private cache copy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "threshold_map"))
import run_map as R5  # noqa: E402
from threebody.coupled import CoupledBasis  # noqa: E402

TOL = 1e-12


class TwoPass(CoupledBasis):
    def project(self, tol=1e-12):
        if self._proj is not None and self._proj["tol"] == tol:
            return self._proj
        pr = CoupledBasis.project(self, tol)
        S = self.mats["S"]
        X = pr["X"]
        lam2, U2 = np.linalg.eigh(X.T @ S @ X)
        keep2 = lam2 > tol * lam2[-1]
        pr.update(X=X @ (U2[:, keep2] / np.sqrt(lam2[keep2])), kept=int(keep2.sum()), P={})
        self._proj = pr
        return pr


def main(targets):
    out = open(os.path.join(HERE, "check_run5_orthog.txt"), "a")
    P = lambda t: (print(t, flush=True), out.write(t + "\n"), out.flush())
    P(f"# check_run5_orthog.py: run-5 standard largest bases, tol {TOL:g}, single vs two-pass orthogonalization")
    for t in targets:
        cfg = R5.config(t)
        sc = json.load(open(os.path.join(R5.DATA, f"{t}_scales.json")))
        mu = sc["M_opt"] / (sc["M_opt"] + 1)
        nd = len(cfg["diffuse"])
        rows = list(csv.DictReader(open(os.path.join(R5.DATA, f"{t}_crossings.csv"))))
        lab = R5.full_label(cfg, cfg["ND"])
        prod = [float(r["M_star"]) for r in rows if r["label"] == lab and float(r["tol"]) == TOL and float(r["margin"]) == 0.0]
        b = R5.make_basis(cfg, sc["scales"], cfg["ND"], nd, mu)
        S = b.mats["S"]
        pr = b.project(TOL)
        X = pr["X"]
        G = X.T @ S @ X
        err1 = float(np.max(np.abs(G - np.eye(G.shape[0]))))
        lam2 = np.linalg.eigvalsh(G)
        lo, hi = R5.bracket(cfg)
        M1 = R5.crossing(cfg, b, lo, hi, TOL)
        b.__class__ = TwoPass
        b._proj = None
        M2 = R5.crossing(cfg, b, lo, hi, TOL)
        P(f"{t} '{lab}' size={b.size} kept={pr['kept']}: single-pass max|X^T S X - I| = {err1:.1e} "
          f"(eig {lam2[0]:.6f}..{lam2[-1]:.6f}); M* single {M1:.10f} (production {prod[0] if prod else float('nan'):.10f}), "
          f"two-pass {M2:.10f}, shift {M2 - M1:+.2e}")
        del b


if __name__ == "__main__":
    main(sys.argv[1:])
