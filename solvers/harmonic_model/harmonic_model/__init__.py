"""Exactly solvable three-body testbeds for the {X, X, x} mass family.

Model A (``model_a``): all-harmonic (Moshinsky type), springs k (unlike pairs)
and k' (like pair).
Model B (``model_b``): harmonic unlike pairs plus g/R^2 repulsion between the
like particles (Calogero type). Ground state in closed form.
``twobody``: two-body harmonic and Coulomb (masses M and 1) reference results.
``numerics``: independent numerical solvers used only for validation.

Units: hbar = 1, mass of the odd particle = 1. Particles 1, 2 have mass M.
Jacobi coordinates R = r2 - r1, rho = r3 - (r1 + r2)/2.

BLAS threads are limited (HARMONIC_BLAS_THREADS, default 2) if numpy has not
been imported yet.
"""
import os as _os

_n = _os.environ.get("HARMONIC_BLAS_THREADS", "2")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, _n)

try:  # enforce the limit even if numpy was loaded first
    from threadpoolctl import threadpool_limits as _tl
    _tl(int(_n))
except Exception:  # pragma: no cover
    pass

from . import common, model_a, model_b, twobody  # noqa: E402,F401

__all__ = ["common", "model_a", "model_b", "twobody"]
