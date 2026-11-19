"""L = 6 / 7 lowest energy vs compact size N and canonical-orthogonalization tol (validate_L.py system: H2+,
M = 1836.15267247, J = L).  Written 2026-11-16 for the L = 7 validation failure; output in
validate_L7_diagnostics.md.  Uses the local term cache read-only.  Run from the project root."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "solvers"))
import tb_local  # noqa  (run 5's cache, read only: all L = 6, 7 lists exist)
from threebody import symmetric_family
from threebody.coupled import CoupledBasis, block
from scipy.linalg import eigh
AL = (9.95811007, 9.95811007, 1.42227562); s = symmetric_family(1836.15267247)
BOAD = {6: -0.591710122507, 7: -0.589964708878}
for L in (6, 7):
    for n in (12, 14, 16, 18, 20):
        t0 = time.time()
        b = CoupledBasis(L, block(n, AL, (1, 1, 2)), "A" if L % 2 else "S")
        for tol in (1e-10, 1e-12, 1e-14):
            H = b.hamiltonian_projected(s, tol)
            e = eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0]
            print(f"L={L} N={n} size={b.size} tol={tol:g} kept={H.shape[0]} E={e:.12f} E-E_BOad={e - BOAD[L]:+.3e} ({time.time()-t0:.0f}s)", flush=True)
