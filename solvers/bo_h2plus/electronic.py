"""Exact clamped-nuclei one-electron two-center problem (homonuclear, charges Z, Z).

Atomic units, electron mass 1 unless stated.  Prolate spheroidal coordinates
xi = (r1 + r2)/R in [1, inf), eta = (r1 - r2)/R in [-1, 1].  For m = 0 the
Schroedinger equation separates into

    d/deta[(1-eta^2) Y'] + (-lam + p^2 eta^2) Y = 0                  (angular)
    d/dxi [(xi^2-1) X'] + (lam + a xi - p^2 xi^2) X = 0              (radial)

with a = 2 Z R and p^2 = -E R^2 / 2, E the electronic energy (without the
nuclear repulsion Z^2/R).

Angular equation: expansion in normalized Legendre polynomials (even l for
sigma_g, odd l for sigma_u); lam is an eigenvalue of a symmetric tridiagonal
matrix (largest one for the lowest state of each parity).

Radial equation: Jaffe expansion
    X = (xi+1)^sigma exp(-p xi) sum_n g_n t^n,  t = (xi-1)/(xi+1),
    sigma = a/(2p) - 1,
which gives the three-term recurrence (derived with sympy, see README)
    (n+1)^2 g_{n+1} + [lam + Bt_n] g_n + (n-1-sigma)^2 g_{n-1} = 0,
    Bt_n = -(2n^2 + 4np + 2n + p^2 + 2p + 1 - (a/p) n - a - a/(2p)).
Truncated, this is a tridiagonal eigenproblem for -lam that can be
symmetrized (off-diagonal (n+1)|n-sigma|); the nodeless radial function is the
one with the smallest lam.

The energy follows from the matching condition lam_rad(p) = lam_ang(p),
solved for p by Brent's method.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq

__all__ = [
    "electronic_energy", "total_energy", "dEdR", "equilibrium",
    "electronic_energy_mu", "Wavefunction", "wavefunction",
]


def _sizes(p, n_ang=None, n_rad=None):
    if n_ang is None:
        n_ang = 30 + int(1.2 * p)
    if n_rad is None:
        n_rad = 60
    return n_ang, n_rad


def _ang_matrix(p, parity, n):
    l = np.arange(parity, parity + 2 * n, 2).astype(float)
    d = -l * (l + 1) + p * p * (2 * l * l + 2 * l - 1) / ((2 * l - 1) * (2 * l + 3))
    ll = l[:-1]
    e = p * p * (ll + 1) * (ll + 2) / ((2 * ll + 3) * np.sqrt((2 * ll + 1) * (2 * ll + 5)))
    return l, d, e


def _lam_ang(p, parity, n, vectors=False):
    l, d, e = _ang_matrix(p, parity, n)
    k = n - 1
    if vectors:
        w, v = eigh_tridiagonal(d, e, select="i", select_range=(k, k))
        return w[0], l, v[:, 0]
    w = eigh_tridiagonal(d, e, eigvals_only=True, select="i", select_range=(k, k))
    return w[0]


def _rad_matrix(p, a, n):
    s = a / (2 * p) - 1.0
    nn = np.arange(n, dtype=float)
    bt = -(2 * nn ** 2 + 4 * nn * p + 2 * nn + p ** 2 + 2 * p + 1 - (a / p) * nn - a - a / (2 * p))
    e = (nn[:-1] + 1) * np.abs(nn[:-1] - s)
    return s, bt, e


def _lam_rad(p, a, n, vectors=False):
    """Smallest lam of the radial problem (= -largest eigenvalue of the matrix)."""
    s, bt, e = _rad_matrix(p, a, n)
    k = n - 1
    if not vectors:
        w = eigh_tridiagonal(bt, e, eigvals_only=True, select="i", select_range=(k, k))
        return -w[0]
    w, v = eigh_tridiagonal(bt, e, select="i", select_range=(k, k))
    v = v[:, 0]
    # undo the symmetrization: g = D v, d_{n+1}/d_n = |n - s|/(n+1)
    nn = np.arange(n - 1, dtype=float)
    ratio = np.abs(nn - s) / (nn + 1)
    logd = np.concatenate([[0.0], np.cumsum(np.log(np.maximum(ratio, 1e-300)))])
    g = v * np.exp(logd)
    return -w[0], s, g


def _match(p, R, Z, parity, n_ang, n_rad):
    na, nr = _sizes(p, n_ang, n_rad)
    return _lam_rad(p, 2 * Z * R, nr) - _lam_ang(p, parity, na)


def _solve_p(R, Z, parity, n_ang, n_rad):
    if R <= 0:
        raise ValueError("R must be > 0")
    # electronic energy lies in (-2 Z^2, -Z^2/8) for both 1s sigma_g and 2p sigma_u
    e_lo, e_hi = -2.0 * Z * Z * (1 + 1e-9), -Z * Z / 8.0
    p_lo, p_hi = R * math.sqrt(-e_hi / 2.0), R * math.sqrt(-e_lo / 2.0)
    f = lambda p: _match(p, R, Z, parity, n_ang, n_rad)
    flo, fhi = f(p_lo), f(p_hi)
    if flo * fhi > 0:
        # scan for a sign change (should not happen for the two supported states)
        ps = np.linspace(p_lo, p_hi, 200)
        fs = [f(x) for x in ps]
        for i in range(len(ps) - 1):
            if fs[i] * fs[i + 1] <= 0:
                p_lo, p_hi = ps[i], ps[i + 1]
                break
        else:
            raise RuntimeError(f"no root bracketed at R={R}")
    return brentq(f, p_lo, p_hi, xtol=1e-16, rtol=1e-15, maxiter=200)


_STATES = {"1sg": 0, "2pu": 1}


def electronic_energy(R, state="1sg", Z=1.0, n_ang=None, n_rad=None):
    """Electronic energy (no nuclear repulsion), electron mass 1, nuclear charges Z.

    state: '1sg' (1s sigma_g, ground) or '2pu' (2p sigma_u).
    """
    parity = _STATES[state]
    R = float(R)
    p = _solve_p(R, Z, parity, n_ang, n_rad)
    return -2.0 * p * p / (R * R)


def total_energy(R, state="1sg", Z=1.0, **kw):
    """BO potential U(R) = E_el(R) + Z^2/R, electron mass 1."""
    return electronic_energy(R, state, Z, **kw) + Z * Z / float(R)


def electronic_energy_mu(R, mu, state="1sg"):
    """Electronic energy for electron (reduced) mass mu: mu * E_el(mu R; 1)."""
    return mu * electronic_energy(mu * R, state)


def dEdR(R, state="1sg", Z=1.0, h=None, total=True):
    """dU/dR (total=True) or dE_el/dR by a Richardson-extrapolated central difference."""
    R = float(R)
    if h is None:
        h = min(2e-3 * max(R, 0.05), 0.05)
    f = (lambda r: total_energy(r, state, Z)) if total else (lambda r: electronic_energy(r, state, Z))
    d1 = (f(R + h) - f(R - h)) / (2 * h)
    d2 = (f(R + h / 2) - f(R - h / 2)) / h
    return (4 * d2 - d1) / 3


@dataclass
class Equilibrium:
    Re: float
    Emin: float
    k: float  # d2U/dR2 at Re (hartree/bohr^2)


def equilibrium(state="1sg", mu=1.0, bracket=(1.5, 2.5)):
    """Minimum of U(R) = mu E_el(mu R) + 1/R; returns R_e, E_min, force constant k.

    Uses the exact mass scaling U_mu(R) = mu U_1(mu R): R_e = R_e1/mu,
    E_min = mu E_min1, k = mu^3 k1.
    """
    if state != "1sg":
        raise NotImplementedError("equilibrium only for 1s sigma_g")
    Re1 = brentq(lambda r: dEdR(r), *bracket, xtol=1e-13)
    E1 = total_energy(Re1)
    h = 1e-2
    f = total_energy
    k1a = (f(Re1 + h) - 2 * E1 + f(Re1 - h)) / h ** 2
    k1b = (f(Re1 + h / 2) - 2 * E1 + f(Re1 - h / 2)) / (h / 2) ** 2
    k1 = (4 * k1b - k1a) / 3
    return Equilibrium(Re1 / mu, mu * E1, mu ** 3 * k1)


class Wavefunction:
    """Normalized 1s sigma_g / 2p sigma_u wavefunction psi(xi, eta) (unit mass, charge Z).

    Normalization: 2 pi (R/2)^3 int int psi^2 (xi^2 - eta^2) dxi deta = 1.
    """

    def __init__(self, R, state="1sg", Z=1.0, n_ang=None, n_rad=None):
        self.R, self.Z, self.state = float(R), Z, state
        parity = _STATES[state]
        self.p = p = _solve_p(self.R, Z, parity, n_ang, n_rad)
        self.E = -2 * p * p / self.R ** 2
        na, nr = _sizes(p, n_ang, n_rad)
        lam_a, self.l, c = _lam_ang(p, parity, na, vectors=True)
        lam_r, self.sigma, g = _lam_rad(p, 2 * Z * self.R, nr, vectors=True)
        self.lam = 0.5 * (lam_a + lam_r)
        # fix signs: Y(1) > 0, X(1) = g_0 > 0
        if np.sum(c * np.sqrt(self.l + 0.5)) < 0:
            c = -c
        if g[0] < 0:
            g = -g
        self.c = c * np.sqrt(self.l + 0.5)  # coefficients of plain P_l
        self.g = g
        self.N = 1.0
        self.N = 1.0 / math.sqrt(self._norm2())

    # --- separated factors -------------------------------------------------
    def X(self, xi, deriv=False):
        xi = np.asarray(xi, float)
        t = (xi - 1) / (xi + 1)
        poly = np.polynomial.polynomial.polyval(t, self.g)
        pref = np.exp(self.sigma * np.log(xi + 1) - self.p * xi)
        if not deriv:
            return pref * poly
        dpoly = np.polynomial.polynomial.polyval(t, np.polynomial.polynomial.polyder(self.g))
        dt = 2.0 / (xi + 1) ** 2
        val = pref * poly
        dval = pref * (poly * (self.sigma / (xi + 1) - self.p) + dpoly * dt)
        return val, dval

    def Y(self, eta, deriv=False):
        eta = np.asarray(eta, float)
        from numpy.polynomial import legendre as L
        coef = np.zeros(int(self.l[-1]) + 1)
        coef[self.l.astype(int)] = self.c
        val = L.legval(eta, coef)
        if not deriv:
            return val
        return val, L.legval(eta, L.legder(coef))

    def __call__(self, xi, eta):
        return self.N * self.X(xi) * self.Y(eta)

    def _norm2(self, nx=120, ne=120):
        xi, wx, eta, we = quad_nodes(self.p, nx, ne)
        X = self.X(xi)
        Y = self.Y(eta)
        XX, EE = np.meshgrid(xi, eta, indexing="ij")
        w = np.outer(wx, we) * (XX ** 2 - EE ** 2)
        return 2 * math.pi * (self.R / 2) ** 3 * np.sum(w * np.outer(X, Y) ** 2)


def quad_nodes(p, nx=120, ne=120):
    """Quadrature in xi (Gauss-Laguerre in 2p(xi-1)) and eta (Gauss-Legendre)."""
    x, wx = np.polynomial.laguerre.laggauss(nx)
    xi = 1 + x / (2 * p)
    wx = wx * np.exp(x) / (2 * p)
    eta, we = np.polynomial.legendre.leggauss(ne)
    return xi, wx, eta, we


def wavefunction(R, state="1sg", Z=1.0):
    return Wavefunction(R, state, Z)
