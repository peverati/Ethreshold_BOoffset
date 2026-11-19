"""Nonadiabatic vibrational and rotational mass functions m_v(R), m_r(R) for {X+, X+, e-}.

Second-order nonadiabatic perturbation theory (see THEORY.md).  Unit electron mass, R in bohr
(this is the scaled coordinate r = mu_e R of the project's lambda = M + 1/2 convention):

    A_par(R)  = < f_par | Q (h - E0)^-1 Q | f_par >,   f_par = d_R phi0 |_rho        (sigma_g -> sigma_g)
    A_perp(R) = < f_perp | (h - E0)^-1 | f_perp >,     f_perp = (1/R) L_y phi0          (sigma_g -> pi_g)
    m_v(R) = 4 A_par(R),   m_r(R) = 4 A_perp(R)        (M -> M + m in lambda),  m(inf) = 1/2.

The resolvent is applied by solving the Sternheimer equation (h - E0) psi1 = -Q f in a
prolate-spheroidal basis (exact Gauss quadrature for all H and S integrals, Hylleraas
functional, i.e. the basis result is a lower bound converging from below):

    m = 0 (sigma_g):  chi_{i l} = L_i(x) e^{-x/2} P_l(eta),                  l = 0, 2, 4, ...
    m = 1 (pi_g)   :  chi_{i l} = sqrt(xi^2-1) L_i(x) e^{-x/2} P_l^1(eta),   l = 2, 4, 6, ...  (times cos phi)
    x = 2 alpha (xi - 1),  alpha = p = R sqrt(-E0/2)  (the exact decay of phi0).

phi0 and its derivatives come from solvers/bo_h2plus (exact to ~1e-12); f_par uses the same
chain rule as bo_h2plus/adiabatic.py.  Basis sizes: n_xi Laguerre functions, n_eta Legendre
functions.  Run nothing here; see run.py and test_massfunc.py.
"""
from __future__ import annotations

import math
import os
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_here, "..", "..", "..", "solvers"))

import numpy as np
from scipy.linalg import eigh

from bo_h2plus.electronic import Wavefunction

try:
    from threadpoolctl import threadpool_limits
except ImportError:  # pragma: no cover
    threadpool_limits = None

NTHREADS = int(os.environ.get("MASSFUNC_THREADS", "4"))


# ---------------------------------------------------------------- 1D basis pieces
def _laguerre(n, x):
    """L_0..L_{n-1}(x) and derivatives, shape (n, len(x)) (stable recurrence)."""
    L = np.zeros((n, x.size))
    L[0] = 1.0
    if n > 1:
        L[1] = 1.0 - x
    for k in range(1, n - 1):
        L[k + 1] = ((2 * k + 1 - x) * L[k] - k * L[k - 1]) / (k + 1)
    dL = np.zeros_like(L)
    dL[1:] = -np.cumsum(L[:-1], axis=0)        # L_n' = -sum_{k<n} L_k
    return L, dL


def _legendre(lmax, t):
    """P_0..P_lmax(t), P', P'' (t strictly inside (-1, 1))."""
    P = np.zeros((lmax + 1, t.size))
    P[0] = 1.0
    if lmax >= 1:
        P[1] = t
    for l in range(1, lmax):
        P[l + 1] = ((2 * l + 1) * t * P[l] - l * P[l - 1]) / (l + 1)
    l = np.arange(lmax + 1)[:, None]
    one = 1.0 - t * t
    dP = np.zeros_like(P)
    dP[1:] = l[1:] * (P[:-1] - t * P[1:]) / one
    d2P = (2 * t * dP - l * (l + 1) * P) / one
    return P, dP, d2P


def _xi_basis(n, alpha, xi, m):
    x = 2 * alpha * (xi - 1)
    L, dL = _laguerre(n, x)
    e = np.exp(-x / 2)
    u = L * e
    du = 2 * alpha * (dL - L / 2) * e
    if m == 0:
        return u, du
    s = np.sqrt(xi * xi - 1)
    return s * u, (xi / s) * u + s * du


def _eta_basis(n, eta, m):
    if m == 0:
        ls = np.arange(0, 2 * n, 2)
    else:
        ls = np.arange(2, 2 * n + 2, 2)
    P, dP, d2P = _legendre(int(ls[-1]), eta)
    P, dP, d2P = P[ls], dP[ls], d2P[ls]
    if m == 0:
        nrm = np.sqrt(ls + 0.5)[:, None]
        return nrm * P, nrm * dP
    s = np.sqrt(1 - eta * eta)
    nrm = np.sqrt((2 * ls + 1) / (2.0 * ls * (ls + 1)))[:, None]
    Y = s * dP
    dY = -(eta / s) * dP + s * d2P
    return nrm * Y, nrm * dY


def _nodes(alpha, nqx, nqe):
    x, wx = np.polynomial.laguerre.laggauss(nqx)
    xi = 1 + x / (2 * alpha)
    wx = np.exp(np.log(wx) + x) / (2 * alpha)           # plain d xi weights
    eta, we = np.polynomial.legendre.leggauss(nqe)
    return xi, wx, eta, we


# ---------------------------------------------------------------- right-hand sides
def _rhs_functions(R, xi, eta, h=None):
    """f_par = d_R phi0|_rho and g = (eta phi_xi - xi phi_eta) sqrt((xi^2-1)(1-eta^2))/(xi^2-eta^2)
    (the rotation generator; L_y phi0 = +-i cos(phi) g), on the (xi, eta) grid.  Unit mass."""
    if h is None:
        h = 1e-3 * R
    wf = Wavefunction(R)
    X, dX = wf.X(xi, deriv=True)
    Y, dY = wf.Y(eta, deriv=True)
    psi = wf.N * np.outer(X, Y)
    pxi = wf.N * np.outer(dX, Y)
    peta = wf.N * np.outer(X, dY)
    XI, ET = np.meshgrid(xi, eta, indexing="ij")
    den = XI ** 2 - ET ** 2
    vals = {}
    for s in (-2, -1, 1, 2):
        w = Wavefunction(R + s * h)
        vals[s] = w.N * np.outer(w.X(xi), w.Y(eta))
    dR_fixed = (8 * (vals[1] - vals[-1]) - (vals[2] - vals[-2])) / (12 * h)
    F = dR_fixed - (XI * (XI ** 2 - 1) * pxi + ET * (1 - ET ** 2) * peta) / (R * den)
    G = np.sqrt((XI ** 2 - 1) * (1 - ET ** 2)) * (ET * pxi - XI * peta) / den
    return wf, psi, F, G


# ---------------------------------------------------------------- one symmetry block
def _block(R, m, n_xi, n_eta, alpha, xi, wx, eta, we):
    U, dU = _xi_basis(n_xi, alpha, xi, m)
    Y, dY = _eta_basis(n_eta, eta, m)
    x2 = xi * xi
    e2 = eta * eta
    mx = lambda wt, A, B: (A * (wx * wt)) @ B.T
    Xa, Xb, Xk, Xv = mx(x2, U, U), mx(1.0, U, U), mx(x2 - 1, dU, dU), mx(xi, U, U)
    Ea, Eb, Ek = (Y * we) @ Y.T, (Y * (we * e2)) @ Y.T, (dY * (we * (1 - e2))) @ dY.T
    S = (R / 2) ** 3 * (np.kron(Xa, Ea) - np.kron(Xb, Eb))
    T = (R / 4) * (np.kron(Xk, Ea) + np.kron(Xb, Ek))
    if m:
        Xc = mx(1.0 / (x2 - 1), U, U)
        Ec = (Y * (we / (1 - e2))) @ Y.T
        T += (R / 4) * m * m * (np.kron(Xb, Ec) + np.kron(Xc, Ea))
    V = -(R * R / 2) * np.kron(Xv, Ea)
    return U, Y, S, T + V


def _project(U, Y, R, wx, we, xi, eta, f):
    """b_k = int chi_k f (R/2)^3 (xi^2-eta^2) dxi deta  (phi integral not included)."""
    XI, ET = np.meshgrid(xi, eta, indexing="ij")
    W = (R / 2) ** 3 * np.outer(wx, we) * (XI ** 2 - ET ** 2)
    return np.einsum("iq,lr,qr->il", U, Y, W * f).ravel()


def default_sizes(R, level=0):
    """Basis sizes (n_xi, n_eta) for a given R. level -1 small, 0 standard, 1 large, 2 very large."""
    p = R * math.sqrt(0.25 + 0.5 / max(R, 0.1))           # rough p ~ R sqrt(-E/2)
    n_xi = {-1: 18, 0: 24, 1: 32, 2: 40}[level]
    n_eta = int(round({-1: 10, 0: 14, 1: 20, 2: 26}[level] + {-1: 1.0, 0: 1.3, 1: 1.8, 2: 2.2}[level] * p))
    return n_xi, n_eta


def mass_functions(R, n_xi=None, n_eta=None, level=0, alpha_scale=1.0, details=True):
    """Return dict with m_v, m_r (and diagnostics) at internuclear distance R (unit mass)."""
    t0 = time.time()
    ctx = threadpool_limits(NTHREADS) if threadpool_limits else None
    try:
        a, b = default_sizes(R, level)
        n_xi = n_xi or a
        n_eta = n_eta or b
        wf0 = Wavefunction(R)
        alpha = alpha_scale * wf0.p
        nqx = n_xi + 40
        nqe = 2 * n_eta + 40
        xi, wx, eta, we = _nodes(alpha, nqx, nqe)
        wf, psi, F, G = _rhs_functions(R, xi, eta)
        E0 = wf.E
        out = dict(R=R, n_xi=n_xi, n_eta=n_eta, alpha=alpha, E0_exact=E0)
        # ---- sigma_g: radial (vibrational) mass
        U, Y, S, H = _block(R, 0, n_xi, n_eta, alpha, xi, wx, eta, we)
        ev, C = eigh(H, S)
        bF = _project(U, Y, R, wx, we, xi, eta, F)
        bP = _project(U, Y, R, wx, we, xi, eta, psi)
        cF = C.T @ bF
        A_par = 2 * math.pi * np.sum(cF[1:] ** 2 / (ev[1:] - ev[0]))
        out.update(E0_basis=ev[0], E1_sigma_g=ev[1],
                   # completeness of the basis for f_par: sum_n <n|f>^2 / <f|f>_radial-part
                   fpar_norm_in_basis=2 * math.pi * np.sum(cF ** 2),
                   fpar_overlap_phi0=2 * math.pi * cF[0],
                   phi0_overlap=2 * math.pi * (C[:, 0] @ bP) ** 2,
                   A_par=A_par, m_v=4 * A_par)
        # ---- pi_g: perpendicular (rotational) mass
        U1, Y1, S1, H1 = _block(R, 1, n_xi, n_eta, alpha, xi, wx, eta, we)
        ev1, C1 = eigh(H1, S1)
        bG = _project(U1, Y1, R, wx, we, xi, eta, G / R)
        cG = C1.T @ bG
        A_perp = math.pi * np.sum(cG ** 2 / (ev1 - E0))
        out.update(E_pi_g=ev1[0], fperp_norm_in_basis=math.pi * np.sum(cG ** 2),
                   A_perp=A_perp, m_r=4 * A_perp)
        out["seconds"] = time.time() - t0
        return out
    finally:
        if ctx is not None:
            ctx.__exit__(None, None, None)


def reference_norms(R):
    """<f_par|f_par> and <f_perp|f_perp> from bo_h2plus.adiabatic (D_radial, D_angular/2)."""
    from bo_h2plus.adiabatic import D_of_R
    d = D_of_R(R, details=True)
    return d["radial"], d["angular"] / 2


if __name__ == "__main__":
    for R in map(float, sys.argv[1:] or ["2.0"]):
        r = mass_functions(R)
        print({k: (round(v, 10) if isinstance(v, float) else v) for k, v in r.items()})
