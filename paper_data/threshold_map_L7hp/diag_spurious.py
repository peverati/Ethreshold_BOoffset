"""Diagnose the spurious bound level of the J7v0 check basis 'compact N=8 + 2 diffuse' (written 2026-11-17).

That basis is a subset of 'compact N=10 + 2 diffuse' (crossing 87.52147) yet showed E_0 - E_th = -1.5e-5 at
M = 86.25 (no bracket).  For both bases: gap E_0 - E_th at M = 86.2486, 87.0, 88.0486 for projection tol 1e-12 ..
1e-15 (single canonical orthogonalization, as in production) and with a two-pass orthogonalization (X1 from S at tol;
S2 = X1^T S X1 re-diagonalized and thresholded again at tol; X = X1 X2), plus lam_min/lam_max and kept counts.
Usage (this folder): HP_NPROC=20 python diag_spurious.py"""
import json
import os

import numpy as np
from scipy.linalg import eigh

import run_map_hp as R
from threebody import family_threshold, symmetric_family

MS = (86.2485555653976, 87.0, 88.0485555653976)
cfg = R.config("J7v0")
sc = json.load(open(os.path.join(R.DATA, "J7v0_scales.json")))
mu = sc["M_opt"] / (sc["M_opt"] + 1)          # as in run_map.run_crossing


def two_pass(S, tol):
    dd = 1.0 / np.sqrt(np.diag(S))
    lam, U = np.linalg.eigh(S * np.outer(dd, dd))
    keep = lam > tol * lam[-1]
    X = (dd[:, None] * U[:, keep]) / np.sqrt(lam[keep])
    S2 = X.T @ S @ X
    lam2, U2 = np.linalg.eigh(S2)
    keep2 = lam2 > tol * lam2[-1]
    return X @ (U2[:, keep2] / np.sqrt(lam2[keep2])), lam2


for N in (8, 10):
    b = R.make_basis(cfg, sc["scales"], N, 2, mu)
    S = b.mats["S"]
    for tol in (1e-12, 1e-13, 1e-14, 1e-15):
        b._proj = None
        pr = b.project(tol)
        g = []
        for M in MS:
            H = b.hamiltonian_projected(symmetric_family(M), tol)
            g.append(eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0] - family_threshold(M))
        X2, lam2 = two_pass(S, tol)
        g2 = []
        for M in MS:
            H2 = X2.T @ b.hamiltonian(symmetric_family(M)) @ X2
            g2.append(eigh(H2, eigvals_only=True, subset_by_index=[0, 0])[0] - family_threshold(M))
        print(f"N={N}+2 diffuse size={b.size} tol={tol:g}: kept={pr['kept']} lam_min/lam_max={pr['lam_min'] / pr['lam_max']:+.1e} "
              f"gap(M)={', '.join(f'{x:+.3e}' for x in g)} | two-pass kept={X2.shape[1]} "
              f"(2nd-pass eig range {lam2[0]:.3e}..{lam2[-1]:.3e}) gap={', '.join(f'{x:+.3e}' for x in g2)}", flush=True)
    del b
