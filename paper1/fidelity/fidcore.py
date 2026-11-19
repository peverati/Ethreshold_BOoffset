"""Fidelity susceptibilities of the exact {X+,X+,e-} ground state (L = 0) versus ln M.

Convention (same as solvers/harmonic_model README sec. 4):
    F(M, M e^h) = |<psi(M)|psi(M e^h)>| = 1 - chi h^2 / 2 + ...,
    chi = <dpsi|dpsi> - |<psi|dpsi>|^2,  dpsi = d psi / d ln M,
overlaps in FIXED Jacobi coordinates R = x1 - x0, rho = x2 - (x0+x1)/2 (mass independent,
so identical to fixed interparticle distances).

* chi_F   : plain (fixed-coordinate) susceptibility.
* chi_inv : maximize over a common dilation of all internal coordinates,
            chi_inv = chi_F - <dpsi|D psi>^2 / ||D psi||^2,  D = x.grad + 3 (6 internal dims).
* chi_sep : maximize over independent dilations of R and rho,
            chi_sep = chi_F - v^T G^-1 v,  v = (<dpsi|D psi>, <dpsi|D_R psi>),
            G = Gram matrix of (D psi, D_R psi)  [span{D psi, D_R psi} = span{D_R psi, D_rho psi}].

Perimetric forms (s_k = r_kj + r_kl - r_jl, r_ij = (s_i+s_j)/2, measure w = r01 r02 r12 / 4):
    D   = s0 d0 + s1 d1 + s2 d2 + 3                 (Euler operator: all r_ij scale together)
    D_R = R.grad_R|_rho + 3/2
        = r01 d_r01 + (r01^2+r02^2-r12^2)/(4 r02) d_r02 + (r01^2-r02^2+r12^2)/(4 r12) d_r12 + 3/2,
    d_r01 = d0+d1-d2, d_r02 = d0+d2-d1, d_r12 = d1+d2-d0.
Every integral is polynomial x exp in s (exact product Gauss-Laguerre) except the two squared
singular pieces of ||D_R psi||^2, done in polar coordinates (s_i, s2) = (t u, t (1-u)):
exact Gauss-Laguerre in t for each u, graded composite Gauss-Legendre in u.

dpsi in a FIXED basis (alphas of the mass_scans run at this M) from the bordered linear system
    [[H - E S, S c], [c^T S, 0]] [x, l] = [-(H' - E' S) c, 0],   H' = dH/dlnM = -(T0+T1)/M.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "3")

import numpy as np  # noqa: E402
from scipy.linalg import eigh, lu_factor, lu_solve  # noqa: E402
from scipy.special import roots_laguerre, roots_legendre  # noqa: E402

from threebody import laguerre  # noqa: E402
from threebody.basis import PerimetricBasis  # noqa: E402
from threebody.operators import overlap_terms  # noqa: E402
from threebody.solver import MatrixSet  # noqa: E402

SCANDATA = os.path.join(ROOT, "paper1", "mass_scans", "data")
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)


# ----------------------------------------------------------------------------- basis
def scan_basis(M):
    """(N, alphas, weights) used by paper1/mass_scans at mass M (largest basis)."""
    fn = os.path.join(SCANDATA, f"L0_M{M:g}.npz")
    if not os.path.exists(fn):
        fn = os.path.join(SCANDATA, f"L0_M{round(M, 3):g}.npz")
    d = np.load(fn)
    b = json.loads(str(d["meta"]))["basis"]
    return b["N"], tuple(b["alphas"]), tuple(b["weights"]), float(d["energies"][0])


class Problem:
    """Matrices of one fixed basis; sub-bases N' <= N are sub-blocks (nested)."""

    def __init__(self, N, alphas, weights=(1, 1, 1)):
        self.basis = PerimetricBasis(N, alphas, "S", weights)
        self.ms = MatrixSet(self.basis)
        self.subs = {}

    def sub(self, Nsub):
        if Nsub not in self.subs:
            b = self.basis
            if Nsub == b.N:
                self.subs[Nsub] = (b, np.arange(b.size))
            else:
                w = np.array(b.weights)
                mask = (b.index @ w) <= Nsub
                bs = PerimetricBasis(Nsub, b.alphas, "S", b.weights)
                assert np.array_equal(bs.index, b.index[mask])
                self.subs[Nsub] = (bs, np.where(mask)[0])
        return self.subs[Nsub]

    def mats(self, Nsub):
        bs, idx = self.sub(Nsub)
        ix = np.ix_(idx, idx)
        S = self.ms.S[ix]
        T01 = (self.ms.T[0] + self.ms.T[1])[ix]
        T2 = self.ms.T[2][ix]
        V = (self.ms.V[(0, 1)] - self.ms.V[(0, 2)] - self.ms.V[(1, 2)])[ix]
        return bs, S, T01, T2, V


def ground(S, T01, T2, V, M, kin_scale=1.0):
    """Lowest eigenpair of H = kin_scale*(T01/M + T2) + V (c^T S c = 1, sign <F_000|psi> > 0)."""
    H = kin_scale * (T01 / M + T2) + V
    d = 1.0 / np.sqrt(np.diag(S))
    e, U = eigh(H * np.outer(d, d), S * np.outer(d, d), subset_by_index=[0, 0])
    c = d * U[:, 0]
    c /= np.sqrt(c @ S @ c)
    if (S @ c)[0] < 0:
        c = -c
    return e[0], c, H


def derivative(S, H, E, c, dH):
    """x = d c / d lambda for dH = dH/dlambda at fixed basis, with c^T S x = 0."""
    n = len(c)
    Sc = S @ c
    dE = c @ dH @ c
    A = np.zeros((n + 1, n + 1))
    A[:n, :n] = H - E * S
    A[:n, n] = Sc
    A[n, :n] = Sc
    rhs = np.zeros(n + 1)
    rhs[:n] = -(dH @ c - dE * Sc)
    d = np.ones(n + 1)
    d[:n] = 1.0 / np.sqrt(np.diag(S))
    sol = np.linalg.solve(A * np.outer(d, d), rhs * d) * d
    return sol[:n], dE


def chi_fd(S, T01, T2, V, M, c0, hs=(0.02, 0.01)):
    """Fixed-basis finite-difference chi_F from overlaps at M e^{+-h}; Richardson over hs."""
    vals = []
    for h in hs:
        cs = []
        for sgn in (+1, -1):
            _, c, _ = ground(S, T01, T2, V, M * np.exp(sgn * h))
            cs.append(abs(c0 @ S @ c))
        F = 0.5 * (cs[0] + cs[1])
        vals.append(2.0 * (1.0 - F) / h ** 2)
    rich = (4.0 * vals[1] - vals[0]) / 3.0 if len(hs) == 2 and hs[0] == 2 * hs[1] else np.nan
    return vals, rich


# ------------------------------------------------------------------- quadrature engine
def _tables(nmax, alpha, s):
    F = laguerre.evaluate(nmax, alpha, s)
    D = laguerre.derivative_matrix(nmax + 1)
    return F, alpha * (D.T @ F)


def _gl_nodes(alpha, q):
    x, w = roots_laguerre(q)
    return x / alpha, w * np.exp(x) / alpha


def _graded_u(npan=14, nper=16, umin=1e-7):
    """Composite Gauss-Legendre on [0,1], geometrically graded towards both ends."""
    half = np.concatenate([[0.0], np.geomspace(umin, 0.5, npan)])
    edges = np.concatenate([half, 1.0 - half[::-1][1:]])
    x, w = roots_legendre(nper)
    U, W = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        U.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo))
        W.append(0.5 * (hi - lo) * w)
    return np.concatenate(U), np.concatenate(W)


def product_grid_values(Cs, basis, extra=8):
    """Product Gauss-Laguerre grid; returns s-coordinates (broadcast), weights, values, gradients."""
    nm, al = basis.nmax, basis.alphas
    S_, W_, F_, dF_ = [], [], [], []
    for m in range(3):
        s, w = _gl_nodes(al[m], nm[m] + extra)
        F, dF = _tables(nm[m], al[m], s)
        S_.append(s); W_.append(w); F_.append(F); dF_.append(dF)
    s0 = S_[0][:, None, None]; s1 = S_[1][None, :, None]; s2 = S_[2][None, None, :]
    W = W_[0][:, None, None] * W_[1][None, :, None] * W_[2][None, None, :]

    def ev(C, a, b, c):
        return np.einsum("abc,ai,bj,ck->ijk", C, a, b, c, optimize=True)

    out = []
    for k, C in enumerate(Cs):
        v = ev(C, F_[0], F_[1], F_[2])
        if k == 0:  # gradients only for psi
            g = (ev(C, dF_[0], F_[1], F_[2]), ev(C, F_[0], dF_[1], F_[2]), ev(C, F_[0], F_[1], dF_[2]))
        else:
            g = None
        out.append((v, g))
    return (s0, s1, s2), W, out


def polar_singular(C, basis, pair, extra=10, npan=14, nper=16):
    """int over the pair-(i,2) polar grid of the singular piece of ||D_R psi||^2.

    pair = 0: int (r01 r12 / (4 r02)) Bt^2 ds,  Bt = (r01^2+r02^2-r12^2)/4 * d_r02 psi
    pair = 1: int (r01 r02 / (4 r12)) Ct^2 ds,  Ct = (r01^2-r02^2+r12^2)/4 * d_r12 psi
    With (s_i, s2) = (t u, t (1-u)), ds_i ds2 = t dt du and r_i2 = t/2, the integrand becomes
    r01 r_j2 Bt^2 / 2, a polynomial in (t, u, s_j) times exp(-lambda(u) t - alpha_j s_j).
    """
    i, j = (0, 1) if pair == 0 else (1, 0)
    nm, al = basis.nmax, basis.alphas
    sj, wj = _gl_nodes(al[j], nm[j] + extra)
    Fj, dFj = _tables(nm[j], al[j], sj)
    u, wu = _graded_u(npan, nper)
    qt = nm[i] + nm[2] + extra
    xt, wt = roots_laguerre(qt)
    # C axes -> (i, 2, j)
    Ct = np.transpose(C, (i, 2, j))
    Tj = np.einsum("acb,bq->acq", Ct, Fj, optimize=True)       # (n_i, n_2, Q_j)
    Tjd = np.einsum("acb,bq->acq", Ct, dFj, optimize=True)
    total = 0.0
    chunk = max(1, 4000 // qt)
    for k0 in range(0, len(u), chunk):
        uu = u[k0:k0 + chunk]; ww = wu[k0:k0 + chunk]
        lam = al[i] * uu + al[2] * (1.0 - uu)
        t = xt[None, :] / lam[:, None]                          # (U, T)
        W = (ww[:, None] * wt[None, :] * np.exp(xt)[None, :] / lam[:, None]).ravel()
        si = (t * uu[:, None]).ravel(); s2 = (t * (1.0 - uu[:, None])).ravel()
        Fi, dFi = _tables(nm[i], al[i], si)
        F2, dF2 = _tables(nm[2], al[2], s2)
        # values at (p, q): psi, d_i, d_2, d_j
        G = lambda A, B, T: np.einsum("ap,cp,acq->pq", A, B, T, optimize=True)  # noqa: E731
        gi = G(dFi, F2, Tj); g2 = G(Fi, dF2, Tj); gj = G(Fi, F2, Tjd)
        r_i2 = 0.5 * (si + s2)[:, None]
        r_ij = 0.5 * (si[:, None] + sj[None, :])
        r_j2 = 0.5 * (sj[None, :] + s2[:, None])
        d_i2 = gi + g2 - gj
        # r01 = r_ij, and the "other" distance r_j2
        Bt = (r_ij ** 2 + r_i2 ** 2 - r_j2 ** 2) / 4.0 * d_i2
        integrand = r_ij * r_j2 * Bt ** 2 / 2.0
        total += np.einsum("p,pq,q->", W, integrand, wj)
    return total


def dilation_quantities(c_tensor, x_tensor, basis, extra=8, polar_kw=None):
    """All scalars needed for chi_inv and chi_sep (unnormalized; divide by n0)."""
    (s0, s1, s2), W, out = product_grid_values([c_tensor, x_tensor], basis, extra)
    (psi, (g0, g1, g2)), (xd, _) = out
    r01 = 0.5 * (s0 + s1); r02 = 0.5 * (s0 + s2); r12 = 0.5 * (s1 + s2)
    w = r01 * r02 * r12 / 4.0
    d01 = g0 + g1 - g2; d02 = g0 + g2 - g1; d12 = g1 + g2 - g0
    Dpsi = s0 * g0 + s1 * g1 + s2 * g2 + 3.0 * psi
    A = r01 * d01 + 1.5 * psi
    Bt = (r01 ** 2 + r02 ** 2 - r12 ** 2) / 4.0 * d02
    Ct = (r01 ** 2 - r02 ** 2 + r12 ** 2) / 4.0 * d12
    I = lambda f: float(np.sum(W * f))  # noqa: E731
    q = {}
    q["n0"] = I(w * psi ** 2)
    q["nx"] = I(w * xd ** 2)
    q["psi_x"] = I(w * psi * xd)
    q["psi_D"] = I(w * psi * Dpsi)
    q["DD"] = I(w * Dpsi ** 2)
    q["pD"] = I(w * xd * Dpsi)
    # <f|D_R psi> = int [w f A + (r01 r12/4) f Bt + (r01 r02/4) f Ct]
    DRdot = lambda f: I(w * f * A + r01 * r12 / 4.0 * f * Bt + r01 * r02 / 4.0 * f * Ct)  # noqa: E731
    q["psi_DR"] = DRdot(psi)
    q["pR"] = DRdot(xd)
    q["DR"] = DRdot(Dpsi)
    q["RR_poly"] = I(w * A ** 2 + 2 * r01 * r12 / 4.0 * A * Bt + 2 * r01 * r02 / 4.0 * A * Ct
                     + 2 * r01 / 4.0 * Bt * Ct)
    kw = polar_kw or {}
    q["RR_sing0"] = polar_singular(c_tensor, basis, 0, **kw)
    q["RR_sing1"] = polar_singular(c_tensor, basis, 1, **kw)
    q["RR"] = q["RR_poly"] + q["RR_sing0"] + q["RR_sing1"]
    return q


def chis(chiF, q):
    n0 = q["n0"]
    DD, RR, DR = q["DD"] / n0, q["RR"] / n0, q["DR"] / n0
    pD, pR = q["pD"] / n0, q["pR"] / n0
    chi_inv = chiF - pD ** 2 / DD
    G = np.array([[DD, DR], [DR, RR]])
    v = np.array([pD, pR])
    chi_sep = chiF - v @ np.linalg.solve(G, v)
    return dict(chi_inv=chi_inv, chi_sep=chi_sep, DD=DD, RR=RR, DR=DR, pD=pD, pR=pR,
                gram_corr=DR / np.sqrt(DD * RR))


# --------------------------------------------------------- exact mixed-scale overlaps
def overlap_tensors(Ca, ala, Cb, alb):
    """<a|b> (unnormalized) for tensors with different alphas (exact, Gauss-Laguerre)."""
    total = 0.0
    cache = {}
    for coef, p, k, l in overlap_terms():
        Ms = []
        for m in range(3):
            key = (m, p[m])
            if key not in cache:
                cache[key] = laguerre.mixed_moment_matrix(Ca.shape[m] - 1, ala[m], Cb.shape[m] - 1, alb[m], p[m])
            Ms.append(cache[key])
        X = np.einsum("abc,ad,be,cf->def", Ca, Ms[0], Ms[1], Ms[2], optimize=True)
        total += coef * np.sum(X * Cb)
    return total


def dilated_overlap(Ca, ala, Cb, alb, sigma):
    """<a|U_sigma b>/(|a||b|), (U_sigma f)(r) = e^{3 sigma} f(e^sigma r): b's alphas times e^sigma."""
    alb_s = tuple(a * np.exp(sigma) for a in alb)
    num = np.exp(3 * sigma) * overlap_tensors(Ca, ala, Cb, alb_s)
    na = overlap_tensors(Ca, ala, Ca, ala)
    nb = overlap_tensors(Cb, alb, Cb, alb)
    return num / np.sqrt(na * nb)
