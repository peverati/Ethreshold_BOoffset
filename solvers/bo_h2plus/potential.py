"""Fast, accurate interpolants for the unit-mass BO potential and adiabatic function.

U(R)  = E_el(R) + 1/R  (1s sigma_g, electron mass 1, hartree)
D(R)  = <grad_R psi|grad_R psi> (see adiabatic.py)

Both are tabulated once on piecewise Chebyshev grids (geometric segments
between R = 0.01 and 200 bohr) and cached on disk in ``_cache/``.  Outside the
table:
  R > 200 : E_el = -1/2 - 1/R - 9/(4R^4) - 15/(2R^6)  (dipole/quadrupole polarization;
            the R^-6 coefficient is checked against the exact solver in the tests)
            D = 1/4 + (D(200) - 1/4) (200/R)^4  (tiny, ~1e-10)
  R < 0.01: E_el = -2 + c (R/0.01)^2 matched at 0.01 (irrelevant there: 1/R dominates)
            D = D(0.01)
For a physical mass M use U_M(R) = mu_e U(mu_e R) and E_ad = mu_e^2 D(mu_e R)/M.
"""
from __future__ import annotations

import os
import numpy as np

from .electronic import electronic_energy

R_LO, R_HI = 0.01, 200.0
NSEG = 15            # 0.01 * 2^15 = 327 > 200: geometric segments ratio (R_HI/R_LO)^(1/NSEG)
NCHEB_U = 40
NCHEB_D = 28
_VERSION = 1
_CACHE = os.path.join(os.path.dirname(__file__), "_cache")


def _edges():
    return R_LO * (R_HI / R_LO) ** (np.arange(NSEG + 1) / NSEG)


def _cheb_nodes(a, b, n):
    k = np.arange(n)
    x = np.cos(np.pi * (k + 0.5) / n)
    return 0.5 * (a + b) + 0.5 * (b - a) * x, x


class PiecewiseCheb:
    def __init__(self, edges, coefs):
        self.edges = np.asarray(edges)
        self.coefs = np.asarray(coefs)

    @classmethod
    def build(cls, f, edges, n):
        coefs = []
        for a, b in zip(edges[:-1], edges[1:]):
            r, x = _cheb_nodes(a, b, n)
            vals = np.array([f(ri) for ri in r])
            coefs.append(np.polynomial.chebyshev.chebfit(x, vals, n - 1))
        return cls(edges, np.array(coefs))

    def __call__(self, r):
        r = np.atleast_1d(np.asarray(r, float))
        idx = np.clip(np.searchsorted(self.edges, r) - 1, 0, len(self.edges) - 2)
        out = np.empty_like(r)
        for s in np.unique(idx):
            m = idx == s
            a, b = self.edges[s], self.edges[s + 1]
            x = (2 * r[m] - a - b) / (b - a)
            out[m] = np.polynomial.chebyshev.chebval(x, self.coefs[s])
        return out


_tables = {}


def _load(name, builder):
    if name in _tables:
        return _tables[name]
    fn = os.path.join(_CACHE, f"{name}_v{_VERSION}.npz")
    if os.path.exists(fn):
        d = np.load(fn)
        t = PiecewiseCheb(d["edges"], d["coefs"])
    else:
        t = builder()
        os.makedirs(_CACHE, exist_ok=True)
        np.savez(fn, edges=t.edges, coefs=t.coefs)
    _tables[name] = t
    return t


def eel_table():
    return _load("eel", lambda: PiecewiseCheb.build(electronic_energy, _edges(), NCHEB_U))


def d_table():
    from .adiabatic import D_of_R
    return _load("dbo", lambda: PiecewiseCheb.build(
        lambda r: D_of_R(r, na=100, nb=100), _edges(), NCHEB_D))


def E_el_interp(R):
    """Unit-mass electronic energy (no 1/R) from the table, with asymptotic extensions."""
    R = np.atleast_1d(np.asarray(R, float))
    out = np.empty_like(R)
    t = eel_table()
    mid = (R >= R_LO) & (R <= R_HI)
    out[mid] = t(R[mid])
    hi = R > R_HI
    rh = R[hi]
    out[hi] = -0.5 - 1 / rh - 9 / (4 * rh ** 4) - 15 / (2 * rh ** 6)
    lo = R < R_LO
    if np.any(lo):
        e0 = t(np.array([R_LO]))[0]
        out[lo] = -2 + (e0 + 2) * (R[lo] / R_LO) ** 2
    return out


def U(R):
    """Unit-mass BO potential E_el(R) + 1/R."""
    R = np.atleast_1d(np.asarray(R, float))
    return E_el_interp(R) + 1 / R


def D(R):
    R = np.atleast_1d(np.asarray(R, float))
    t = d_table()
    out = np.empty_like(R)
    mid = (R >= R_LO) & (R <= R_HI)
    out[mid] = t(R[mid])
    d_hi = t(np.array([R_HI]))[0]
    hi = R > R_HI
    out[hi] = 0.25 + (d_hi - 0.25) * (R_HI / R[hi]) ** 4
    lo = R < R_LO
    out[lo] = t(np.array([R_LO]))[0]
    return out


def mu_e(M):
    return 2.0 * M / (2.0 * M + 1.0)


def V_BO(R, M):
    """BO potential for particle mass M (electron reduced mass mu_e): mu_e U(mu_e R)."""
    m = mu_e(M)
    return m * U(m * np.asarray(R, float))


def E_ad(R, M):
    """Diagonal adiabatic correction (hartree) for particle mass M."""
    m = mu_e(M)
    return m * m * D(m * np.asarray(R, float)) / M
