"""Second-order nonadiabatic mass matrices of the {1s sigma_g, 2p sigma_u} model space for {A+, B+, e-}.

Written 2026-11-16.  Extension of paper_data/bo_offset/massfunc/massfunc.py (imported, not modified) from the
single g channel to the two-channel (g, u) model of THEORY.md section 4.  Unit electron mass, R = scaled r.

At fixed rho (nuclear-CM electron coordinate) the nuclear derivatives of phi_i (i = g, u) are

    f_par_i  = F_i + beta Z_i          F = d_R phi|rho_m, Z = d_z phi
    f_perp_i = cos(phi) (G_i/R + beta S_i)   G = s d_z phi - z d_s phi, S = d_s phi

(fields from solvers/bo_h2plus/unequal.derivative_fields).  Quasi-degenerate (Loewdin) second order with the
model space {g, u} excluded from the intermediate states (Q = 1 - |g><g| - |u><u|):

    A_ij(e) = <f_i| Q (h - e)^-1 Q |f_j>,     A_ij = [A_ij(e_i) + A_ij(e_j)]/2,     m_ij = 4 A_ij.

By parity the pieces are (Sg, Su, Pg, Pu = Sigma_g, Sigma_u, Pi_g, Pi_u blocks):

    par : A_gg = A[F_g,F_g;Sg';e_g] + b^2 A[Z_g,Z_g;Su';e_g]      A_uu = A[F_u,F_u;Su';e_u] + b^2 A[Z_u,Z_u;Sg';e_u]
          A_gu = b * avg_e (A[F_g,Z_u;Sg';e] + A[Z_g,F_u;Su';e])
    perp: same with G/R and S on Pg, Pu (no exclusion: g and u are Sigma states).

(' = lowest root of the block removed: 1s sigma_g in Sg, 2p sigma_u in Su.)  For beta = 0 the g row is
exactly massfunc's m_v (par) and m_r (perp).  Asymptotically (atoms): m_gg = m_uu = 1/2 + 2 beta^2,
m_gu = -2 beta (par and perp), so the B-localized open channel has m = 2/(1+r)^2 (exact reduced-mass value).

Stored per R are the beta-independent coefficients:
    mv_gg0, mv_gg2, mv_uu0, mv_uu2, mv_gu1   (par)   m_gg = mv_gg0 + b^2 mv_gg2, m_gu = b mv_gu1, ...
    mr_gg0, mr_gg2, mr_uu0, mr_uu2, mr_gu1   (perp)
"""
from __future__ import annotations

import math
import os
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
sys.path.insert(0, os.path.join(ROOT, "paper_data", "bo_offset", "massfunc"))

import numpy as np  # noqa: E402
from scipy.linalg import eigh  # noqa: E402

import massfunc as MF  # noqa: E402  (symmetric-family code, read only)
from bo_h2plus.electronic import Wavefunction  # noqa: E402
from bo_h2plus.unequal import derivative_fields  # noqa: E402


def _eta_basis_p(n, eta, m, parity):
    """Normalized Legendre (m = 0) or associated Legendre P_l^1 (m = 1) functions of given inversion parity.

    gerade: l even (m = 0: 0, 2, ...; m = 1: 2, 4, ...); ungerade: l odd (1, 3, ...)."""
    if m == 0:
        ls = np.arange(0, 2 * n, 2) if parity == "g" else np.arange(1, 2 * n + 1, 2)
    else:
        ls = np.arange(2, 2 * n + 2, 2) if parity == "g" else np.arange(1, 2 * n + 1, 2)
    P, dP, d2P = MF._legendre(int(ls[-1]), eta)
    P, dP, d2P = P[ls], dP[ls], d2P[ls]
    if m == 0:
        nrm = np.sqrt(ls + 0.5)[:, None]
        return nrm * P, nrm * dP
    s = np.sqrt(1 - eta * eta)
    nrm = np.sqrt((2 * ls + 1) / (2.0 * ls * (ls + 1)))[:, None]
    return nrm * s * dP, nrm * (-(eta / s) * dP + s * d2P)


def _block(R, m, parity, n_xi, n_eta, alpha, xi, wx, eta, we):
    """Same matrices as massfunc._block, for either inversion parity."""
    U, dU = MF._xi_basis(n_xi, alpha, xi, m)
    Y, dY = _eta_basis_p(n_eta, eta, m, parity)
    x2, e2 = xi * xi, eta * eta
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
    ev, C = eigh(T + V, S)
    return U, Y, ev, C


def mass_matrices(R, level=0, alpha_scale=1.0, n_xi=None, n_eta=None):
    """beta-independent mass-matrix coefficients at R (unit mass).  See module docstring."""
    t0 = time.time()
    ctx = MF.threadpool_limits(MF.NTHREADS) if MF.threadpool_limits else None
    try:
        a, b = MF.default_sizes(R, level)
        n_xi, n_eta = n_xi or a, n_eta or b
        wg0, wu0 = Wavefunction(R, "1sg"), Wavefunction(R, "2pu")
        e_g, e_u = wg0.E, wu0.E
        alpha = alpha_scale * min(wg0.p, wu0.p)
        xi, wx, eta, we = MF._nodes(alpha, n_xi + 40, 2 * n_eta + 40)
        _, g = derivative_fields(R, xi, eta, "1sg")
        _, u = derivative_fields(R, xi, eta, "2pu")
        out = dict(R=R, n_xi=n_xi, n_eta=n_eta, alpha=alpha, e_g=e_g, e_u=e_u)
        proj = lambda U, Y, f: MF._project(U, Y, R, wx, we, xi, eta, f)

        def resolvent(block, f1, f2, e, excl, fac):
            U, Y, ev, C = block
            c1, c2 = C.T @ proj(U, Y, f1), C.T @ proj(U, Y, f2)
            k = 1 if excl else 0
            return fac * float(np.sum(c1[k:] * c2[k:] / (ev[k:] - e)))

        # ---- parallel (Sigma blocks, phi integral 2 pi)
        Sg = _block(R, 0, "g", n_xi, n_eta, alpha, xi, wx, eta, we)
        Su = _block(R, 0, "u", n_xi, n_eta, alpha, xi, wx, eta, we)
        out.update(E_Sg0=Sg[2][0], E_Sg1=Sg[2][1], E_Su0=Su[2][0], E_Su1=Su[2][1])
        tp = 2 * math.pi
        A = lambda blk, f1, f2, e: resolvent(blk, f1, f2, e, True, tp)
        out["Apar_gg0"] = A(Sg, g["F"], g["F"], e_g)
        out["Apar_gg2"] = A(Su, g["Z"], g["Z"], e_g)
        out["Apar_uu0"] = A(Su, u["F"], u["F"], e_u)
        out["Apar_uu2"] = A(Sg, u["Z"], u["Z"], e_u)
        out["Apar_gu1"] = 0.5 * sum(A(Sg, g["F"], u["Z"], e) + A(Su, g["Z"], u["F"], e) for e in (e_g, e_u))
        # ---- perpendicular (Pi blocks, cos^2 phi integral pi)
        Pg = _block(R, 1, "g", n_xi, n_eta, alpha, xi, wx, eta, we)
        Pu = _block(R, 1, "u", n_xi, n_eta, alpha, xi, wx, eta, we)
        out.update(E_Pg0=Pg[2][0], E_Pu0=Pu[2][0])
        Ap = lambda blk, f1, f2, e: resolvent(blk, f1, f2, e, False, math.pi)
        Gg, Gu = g["G"] / R, u["G"] / R
        out["Aperp_gg0"] = Ap(Pg, Gg, Gg, e_g)
        out["Aperp_gg2"] = Ap(Pu, g["S"], g["S"], e_g)
        out["Aperp_uu0"] = Ap(Pu, Gu, Gu, e_u)
        out["Aperp_uu2"] = Ap(Pg, u["S"], u["S"], e_u)
        out["Aperp_gu1"] = 0.5 * sum(Ap(Pg, Gg, u["S"], e) + Ap(Pu, g["S"], Gu, e) for e in (e_g, e_u))
        for k in ("gg0", "gg2", "uu0", "uu2", "gu1"):
            out[f"mv_{k}"] = 4 * out[f"Apar_{k}"]
            out[f"mr_{k}"] = 4 * out[f"Aperp_{k}"]
        out["seconds"] = time.time() - t0
        return out
    finally:
        if ctx is not None:
            ctx.__exit__(None, None, None)


KEYS = ["mv_gg0", "mv_gg2", "mv_uu0", "mv_uu2", "mv_gu1", "mr_gg0", "mr_gg2", "mr_uu0", "mr_uu2", "mr_gu1"]


def default_grid():
    """R grid for the mass-matrix table (unit-mass bohr)."""
    return np.unique(np.concatenate([np.geomspace(0.05, 1.0, 15), np.arange(1.1, 12.0, 0.1),
                                     np.arange(12.25, 25.0, 0.25), np.arange(25.5, 40.01, 0.5)]))


if __name__ == "__main__":
    for R in map(float, sys.argv[1:] or ["2.0"]):
        q = mass_matrices(R)
        print({k: round(q[k], 8) for k in KEYS + ["seconds"]})
