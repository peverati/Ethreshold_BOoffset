"""Step 1 diagnostic: double vs double-double assembly for the validate_L basis (H2+, alphas AL, weights (1,1,2)).
Written 2026-11-17.  Usage: python diag_precision.py L N [N ...]"""
import sys
import time

import numpy as np
from scipy.linalg import eigh

import tb_hp_local  # noqa: F401
from hp_assembly import HPCoupledBasis
from threebody import symmetric_family
from threebody.coupled import CoupledBasis, block

AL = (9.95811007, 9.95811007, 1.42227562)
SYS = symmetric_family(1836.15267247)
BOAD = {6: -0.591710122507, 7: -0.589964708878}
L = int(sys.argv[1])
for N in [int(x) for x in sys.argv[2:]]:
    sym = "A" if L % 2 else "S"
    t0 = time.time(); bd = CoupledBasis(L, block(N, AL, (1, 1, 2)), sym); td = time.time() - t0
    t0 = time.time(); bh = HPCoupledBasis(L, block(N, AL, (1, 1, 2)), sym); th = time.time() - t0
    print(f"L={L} N={N} size={bd.size}: build double {td:.0f}s, dd {th:.0f}s", flush=True)
    for op in bd.mats:
        A, B = bd.mats[op], bh.mats[op]
        print(f"  {op:8s} max|double-dd|/max|dd| = {np.max(np.abs(A - B)) / np.max(np.abs(B)):.1e}", flush=True)
    for name, b in (("double", bd), ("dd", bh)):
        S = b.mats["S"]; dd_ = 1 / np.sqrt(np.diag(S)); lam = np.linalg.eigvalsh(S * np.outer(dd_, dd_))
        row = [f"{name:6s} lam_min/lam_max={lam[0] / lam[-1]:+.1e}"]
        for tol in (1e-10, 1e-12, 1e-14, 1e-15):
            b._proj = None
            H = b.hamiltonian_projected(SYS, tol)
            e = eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0]
            row.append(f"tol {tol:g}: kept {H.shape[0]} E-E_BOad {e - BOAD.get(L, 0):+.3e}")
        print("  " + " | ".join(row), flush=True)
