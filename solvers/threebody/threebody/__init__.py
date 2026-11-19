"""Deterministic variational solver for S states of three-particle Coulomb systems.

Perimetric-coordinate Laguerre (Pekeris-type) basis with exact matrix elements.
"""
from .system import ThreeBodySystem, symmetric_family, helium_like, family_threshold, M_PROTON, M_MUON
from .basis import PerimetricBasis
from .solver import MatrixSet, Result, solve, overlap_between
from .convergence import optimize_alphas, convergence_table, extrapolate, default_alphas

__all__ = ["ThreeBodySystem", "symmetric_family", "helium_like", "family_threshold", "M_PROTON", "M_MUON",
           "PerimetricBasis", "MatrixSet", "Result", "solve", "overlap_between",
           "optimize_alphas", "convergence_table", "extrapolate", "default_alphas"]

# OpenBLAS with all cores is badly oversubscribed for the matrix sizes used
# here (a 1500x1500 eigh took ~5 s with 10 threads vs ~0.4 s with 4).
import os as _os

try:
    from threadpoolctl import threadpool_limits as _tl

    _tl(int(_os.environ.get("THREEBODY_BLAS_THREADS", "4")))
except Exception:  # pragma: no cover
    pass
