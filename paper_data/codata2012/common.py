"""Shared settings and helpers for the exact-versus-BO mass scans of {X+, X+, e-}.

Imports the two validated packages unchanged (solvers/threebody, solvers/bo_h2plus).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "3")
os.environ.setdefault("BO_H2PLUS_THREADS", "3")

import numpy as np  # noqa: E402
from scipy.special import roots_legendre  # noqa: E402

from threebody import laguerre  # noqa: E402

DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)

# CODATA 2012 (same as King, Longford & Cox 2013); mass_scans/common.py has CODATA 2018 (1836.15267343, 206.768283)
M_PROTON = 1836.15267245
M_MUON = 206.7682843
# log-spaced backbone plus the requested special masses; denser in 0.3-3 and near critical masses
MASSES = sorted([M_MUON, M_PROTON])   # physical-mass points only
L1_CRIT = 4.2082827   # exact, paper_data/rotational
L2_CRIT = 11.4894217
V1_CRIT = 9.525634    # exact, paper_data/critical_masses
HARTREE_CM = 219474.6313705

# fixed grids for g(R)
R_FIXED = np.linspace(0.0, 20.0, 401)
X_SCALED = np.linspace(0.0, 4.0, 201)   # R/<R>


def kappa(M):
    return M ** -0.25


def family_threshold(M):
    return -M / (2.0 * (M + 1.0))


# ----------------------------------------------------------------------------------------
# One-particle density of X about the centre of mass (L = 0 wave function)
# ----------------------------------------------------------------------------------------
def _composite_gl(a, b, npan, nper, log=True):
    x, w = roots_legendre(nper)
    edges = np.geomspace(a, b, npan + 1) if log else np.linspace(a, b, npan + 1)
    X, W = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        X.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo))
        W.append(0.5 * (hi - lo) * w)
    return np.concatenate(X), np.concatenate(W)


def psi_at(res, state, r01, r02, r12, C=None):
    """Unnormalized wave function (solver convention) at given interparticle distances."""
    b = res.basis
    if C is None:
        C = res.tensor(state)
    s0 = r01 + r02 - r12
    s1 = r01 + r12 - r02
    s2 = r02 + r12 - r01
    F0 = laguerre.evaluate(b.nmax[0], b.alphas[0], s0)
    F1 = laguerre.evaluate(b.nmax[1], b.alphas[1], s1)
    F2 = laguerre.evaluate(b.nmax[2], b.alphas[2], s2)
    T = np.einsum("abc,cq->abq", C, F2, optimize=True)
    T = np.einsum("abq,bq->aq", T, F1, optimize=True)
    return np.einsum("aq,aq->q", T, F0, optimize=True)


def density_about_cm(res, M, s_values, state=0, npan=24, nper=10, ntheta=40, rmax=None):
    """rho_X(s): spherically averaged density of particle X (index 0) about the centre of mass.

    Exact reduction (derivation in README): with R = x0 - x1, r = x2 - (x0+x1)/2 and
    c = 2/(2M+1), the vector from the CM to X0 is s = R/2 - r/(2M+1).  Changing variables
    (R, r) -> (s, r) (Jacobian 8):
        rho(s) = 8 * int d^3r |Psi(R = 2 s zhat + c r, r)|^2,
    independent of the direction zhat because the state has L = 0.  Here |Psi|^2 is the
    probability density in d^3R d^3r = 8 pi^2 r01 r02 r12 dr01 dr02 dr12 (unit Jacobian
    from the solver's relative vectors), i.e. |psi|^2 / (8 pi^2 * res.norm()).
    The remaining 2D integral (r, theta) is done by composite Gauss-Legendre in log r and
    Gauss-Legendre in theta (theta, not cos theta, makes the endpoint coalescence kinks
    smooth).  Returns rho(s) with int 4 pi s^2 rho ds = 1.
    """
    b = res.basis
    C = res.tensor(state)
    nrm = res.norm(state)
    c = 2.0 / (2.0 * M + 1.0)
    if rmax is None:
        rmax = 80.0 / min(b.alphas)
    r, wr = _composite_gl(rmax * 1e-7, rmax, npan, nper, log=True)
    th, wt = roots_legendre(ntheta)
    th = 0.5 * np.pi * (th + 1.0)
    wt = 0.5 * np.pi * wt * np.sin(th)
    rr, tt = np.meshgrid(r, th, indexing="ij")
    W = (wr[:, None] * r[:, None] ** 2 * wt[None, :]).ravel()
    u = np.cos(tt).ravel()
    rr = rr.ravel()
    out = []
    for s in np.atleast_1d(s_values):
        R2 = 4 * s * s + c * c * rr * rr + 4 * s * c * rr * u
        r01 = np.sqrt(np.maximum(R2, 0.0))
        a = 1.0 - c / 2.0
        bb = 1.0 + c / 2.0
        r02 = np.sqrt(np.maximum(rr * rr * a * a + s * s - 2 * s * rr * a * u, 0.0))
        r12 = np.sqrt(np.maximum(rr * rr * bb * bb + s * s + 2 * s * rr * bb * u, 0.0))
        p = psi_at(res, state, r01, r02, r12, C)
        val = 8.0 * 2.0 * np.pi * np.sum(W * p * p) / (8.0 * np.pi ** 2 * nrm)
        out.append(val)
    return np.array(out)


def bo_imports():
    import bo_h2plus  # noqa: F401
    from bo_h2plus import solve_levels, thresholds
    return solve_levels, thresholds
