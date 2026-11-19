"""Electron-nuclear entanglement (R | rho bipartition) of L = 0 three-body states.

*Last updated: 2026-11-13*

Jacobi coordinates: R = x1 - x0 (X-X vector, "nuclear"), rho = x2 - (x0+x1)/2 (electron about the
X-X midpoint, "electronic").  r01 = R, r02^2 = rho^2 + R^2/4 + rho R u, r12^2 = rho^2 + R^2/4 - rho R u,
u = cos(angle(R, rho)).

An L = 0 state is  Psi(R, rho) = sum_l g_l(R, rho) sum_m Y_lm(R^) Y_lm*(rho^),
g_l(R, rho) = 2 pi int_{-1}^{1} du Psi P_l(u)   (Psi normalized in d3R d3rho).
=> the Schmidt decomposition is block diagonal in l; block l has the radial singular values of
g_l (measure R^2 dR, rho^2 drho), each with multiplicity 2l + 1.  Nyström discretization:
A^l_ij = sqrt(a_i) g_l(R_i, rho_j) sqrt(b_j), a_i = wR_i R_i^2, b_j = wrho_j rho_j^2; lambda = sigma^2.

Exact identity used for the orientational / radial split (w_l = (2l+1) ||A^l||_F^2, p_{k|l} the
normalized spectrum of block l):
    S = S_orient + S_rad,  S_orient = -sum_l w_l ln(w_l/(2l+1)),  S_rad = sum_l w_l S_l.
Radial-only state after the LOCAL partial traces over R^ and rho^:  rho_rad = sum_l (2l+1)|g_l><g_l|
(rank <= number of l); its logarithmic negativity is computed in a compressed local basis.
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "3")
os.environ.setdefault("BO_H2PLUS_THREADS", "3")

import numpy as np  # noqa: E402
from scipy.special import roots_legendre, eval_legendre  # noqa: E402

from threebody import laguerre  # noqa: E402
from threebody.basis import PerimetricBasis  # noqa: E402
from threebody.solver import Result  # noqa: E402

SCANDATA = os.path.join(ROOT, "paper_data", "mass_scans", "data")
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)


# ----------------------------------------------------------------------------- grids
def composite_gl(edges, nper):
    x, w = roots_legendre(nper)
    X, W = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        X.append(0.5 * (hi - lo) * x + 0.5 * (hi + lo))
        W.append(0.5 * (hi - lo) * w)
    return np.concatenate(X), np.concatenate(W)


class Grid:
    """R nodes on [Rlo, Rhi]; rho nodes on [0, rmax] (panels uniform in sqrt(rho)); theta GL."""

    def __init__(self, Rlo, Rhi, rmax, nR=(16, 10), nr=(24, 12), nth=128, symmetric=True, lmax=64):
        self.R, self.wR = composite_gl(np.linspace(Rlo, Rhi, nR[0] + 1), nR[1])
        self.r, self.wr = composite_gl(rmax * np.linspace(0, 1, nr[0] + 1) ** 2, nr[1])
        t, wt = roots_legendre(nth)
        top = 0.5 * np.pi if symmetric else np.pi
        th = 0.5 * top * (t + 1.0)
        self.u = np.cos(th)
        self.wu = 0.5 * top * wt * np.sin(th)          # du weights on u in [0,1] or [-1,1]
        self.symmetric = symmetric
        self.ls = np.arange(0, lmax + 1, 2 if symmetric else 1)
        self.P = np.array([eval_legendre(l, self.u) for l in self.ls])   # (nl, nu)
        self.a = self.wR * self.R ** 2
        self.b = self.wr * self.r ** 2
        self.desc = dict(Rlo=Rlo, Rhi=Rhi, rmax=rmax, nR=list(nR), nr=list(nr), nth=nth, lmax=lmax,
                         symmetric=symmetric)


def distances(R, r, u):
    rr, uu = np.meshgrid(r, u, indexing="ij")
    q = rr * rr + 0.25 * R * R
    r02 = np.sqrt(np.maximum(q + rr * R * uu, 0.0))
    r12 = np.sqrt(np.maximum(q - rr * R * uu, 0.0))
    return np.full_like(r02, R), r02, r12


def partial_waves(amp, grid):
    """amp(R, r, u) -> normalized Psi on the (r, u) mesh at fixed R.
    Returns A (nl, nR, nr) and the direct quadrature norm per R node."""
    nl = len(grid.ls)
    A = np.empty((nl, len(grid.R), len(grid.r)))
    direct = np.empty(len(grid.R))
    fac = 4 * np.pi if grid.symmetric else 2 * np.pi        # u in [0,1] counts twice for even l
    nfac = 16 * np.pi ** 2 if grid.symmetric else 8 * np.pi ** 2
    sa, sb = np.sqrt(grid.a), np.sqrt(grid.b)
    for i, R in enumerate(grid.R):
        psi = amp(R, grid.r, grid.u)                       # (nr, nu)
        g = fac * (psi * grid.wu[None, :]) @ grid.P.T      # (nr, nl)
        A[:, i, :] = (g.T * sb[None, :]) * sa[i]
        direct[i] = grid.a[i] * np.sum(grid.b[:, None] * psi * psi * grid.wu[None, :]) * nfac
    return A, direct


# ----------------------------------------------------------------------------- analysis
def _H(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def spectrum(A, ls, nkeep=40):
    """Singular values per l block; returns dict with all entanglement measures."""
    sig = [np.linalg.svd(A[k], compute_uv=False) for k in range(len(ls))]
    lam = [s * s for s in sig]
    mult = 2 * ls + 1
    w = np.array([m * l.sum() for m, l in zip(mult, lam)])
    tot = w.sum()
    w = w / tot
    lam = [l / tot for l in lam]
    allv = np.concatenate([np.repeat(l, m) for l, m in zip(lam, mult)])
    S = _H(allv)
    purity = float(np.sum(allv ** 2))
    Sl = np.array([_H(l / l.sum()) if l.sum() > 0 else 0.0 for l in lam])
    S_orient = _H(w) + float(np.sum(w * np.log(mult)))
    S_rad = float(np.sum(w * Sl))
    lam_sorted = np.sort(allv)[::-1]
    return dict(S=S, S_lin=1 - purity, purity=purity, S_orient=S_orient, S_rad=S_rad,
                w=w, S_l=Sl, lam_max=float(lam_sorted[0]), lam_sorted=lam_sorted[:nkeep],
                lam_blocks=np.array([np.pad(l[:nkeep], (0, max(0, nkeep - len(l[:nkeep])))) for l in lam]),
                norm_sum=float(tot), K=1 / purity)


def radial_measures(A, ls, maxdim=3000):
    """Log-negativity of rho_rad = sum_l (2l+1)|g_l><g_l| and classical mutual information I(|R|;|rho|).
    rho_rad is projected on local bases (R side: span of the block left singular vectors, rho side:
    right ones) keeping singular values with (2l+1) sigma^2 > cut; cut is raised by 10x until
    dim <= maxdim.  The discarded Schmidt weight is returned as 'discard'."""
    mult = 2 * ls + 1
    tot = sum(m * np.sum(a * a) for m, a in zip(mult, A))
    P = sum(m * a * a for m, a in zip(mult, A)) / tot
    pi, qj = P.sum(1), P.sum(0)
    msk = P > 1e-300
    MI = float(np.sum(P[msk] * np.log(P[msk] / np.outer(pi, qj)[msk])))
    svd = [np.linalg.svd(a / math.sqrt(tot), full_matrices=False) for a in A]
    cut = 1e-13
    while True:
        Us, Vs, disc = [], [], 0.0
        for m, (U, s, Vt) in zip(mult, svd):
            keep = m * s * s > cut
            disc += float(np.sum(m * s[~keep] ** 2))
            Us.append(U[:, keep])
            Vs.append(Vt[keep].T)

        def orth(X):
            X = np.concatenate(X, axis=1)
            Q, s, _ = np.linalg.svd(X, full_matrices=False)
            return Q[:, s > 1e-8]
        U = orth(Us)
        V = orth(Vs)
        k, kp = U.shape[1], V.shape[1]
        if k * kp <= maxdim:
            break
        cut *= 10
    G = [U.T @ (a / math.sqrt(tot)) @ V for a in A]
    rho = np.zeros((k * kp, k * kp))
    for m, g in zip(mult, G):
        v = g.reshape(-1)
        rho += m * np.outer(v, v)
    tr = float(np.trace(rho))
    rT = rho.reshape(k, kp, k, kp).transpose(0, 3, 2, 1).reshape(k * kp, k * kp)
    ev = np.linalg.eigvalsh(rT)
    n1 = float(np.sum(np.abs(ev)))
    return dict(MI=MI, E_N=math.log(n1 / tr), neg=(n1 / tr - 1) / 2, k=k, kp=kp, trace=tr, cut=cut,
                discard=disc)


# ----------------------------------------------------------------------------- exact states
def psi_tensor(C, alphas, r01, r02, r12):
    s0 = r01 + r02 - r12
    s1 = r01 + r12 - r02
    s2 = r02 + r12 - r01
    sh = s0.shape
    F0 = laguerre.evaluate(C.shape[0] - 1, alphas[0], s0.ravel())
    F1 = laguerre.evaluate(C.shape[1] - 1, alphas[1], s1.ravel())
    F2 = laguerre.evaluate(C.shape[2] - 1, alphas[2], s2.ravel())
    T = np.tensordot(C, F2, axes=([2], [0]))           # (a, b, q)
    T = np.einsum("abq,bq->aq", T, F1, optimize=True)
    return np.einsum("aq,aq->q", T, F0, optimize=True).reshape(sh)


def scan_file(M):
    for key in (f"{M:g}", f"{round(M, 3):g}"):
        fn = os.path.join(SCANDATA, f"L0_M{key}.npz")
        if os.path.exists(fn):
            return fn
    raise FileNotFoundError(M)


def load_exact(M, state=0):
    """(Result, C tensor, norm, meta-basis) from the mass_scans npz (largest basis)."""
    d = np.load(scan_file(M))
    b = json.loads(str(d["meta"]))["basis"]
    basis = PerimetricBasis(b["N"], tuple(b["alphas"]), "S", tuple(b["weights"]))
    res = Result(system=None, basis=basis, energies=d["energies"], vectors=d["vectors"])
    return res, b


def exact_amp(res, state=0):
    C = res.tensor(state)
    nrm = res.norm(state)
    al = res.basis.alphas
    f = 1.0 / math.sqrt(8 * np.pi ** 2 * nrm)

    def amp(R, r, u):
        r01, r02, r12 = distances(R, r, u)
        return f * psi_tensor(C, al, r01, r02, r12)
    return amp


def exact_amp_prevN(M, dN=4):
    """Ground state in the nested N - dN sub-basis (same alphas), via paper_data/fidelity/fidcore."""
    sys.path.insert(0, os.path.join(ROOT, "paper_data", "fidelity"))
    import fidcore
    _, b = load_exact(M)
    pb = fidcore.Problem(b["N"], tuple(b["alphas"]), tuple(b["weights"]))
    bs, S, T01, T2, V = pb.mats(b["N"] - dN)
    E, c, _ = fidcore.ground(S, T01, T2, V, M)
    res = Result(system=None, basis=bs, energies=np.array([E]), vectors=c[:, None])
    return exact_amp(res), E


# ----------------------------------------------------------------------------- BO product state
def bo_amp(M):
    """Psi_BO = sqrt(P_BO(R)/(4 pi R^2)) phi_1sg(rho; R) (clamped-nuclei electron, BO nuclear J=0 v=0)."""
    from bo_h2plus.nuclear import solve_levels
    from bo_h2plus.electronic import Wavefunction
    ls = solve_levels(M, J=0, adiabatic=False, check=False)
    Rg, P0 = ls.R, ls.P[0]
    cache = {}

    def amp(R, r, u):
        if R not in cache:
            cache[R] = Wavefunction(R)
        wf = cache[R]
        _, ra, rb = distances(R, r, u)
        xi = (ra + rb) / R
        eta = (ra - rb) / R
        p = max(np.interp(R, Rg, P0, left=0.0, right=0.0), 0.0)
        return math.sqrt(p / (4 * np.pi * R * R)) * wf(np.maximum(xi, 1.0), np.clip(eta, -1, 1))
    return amp, ls


def electronic_l_weights(R, rmax=40.0, nr=(40, 14), nth=200, lmax=80):
    """Partial-wave weights of phi_1sg(rho; R) about the midpoint: the M -> inf orientational limit."""
    from bo_h2plus.electronic import Wavefunction
    wf = Wavefunction(R)
    r, wr = composite_gl(rmax * np.linspace(0, 1, nr[0] + 1) ** 2, nr[1])
    # put a panel edge exactly at the nucleus R/2
    edges = np.unique(np.concatenate([rmax * np.linspace(0, 1, nr[0] + 1) ** 2, [R / 2]]))
    r, wr = composite_gl(edges, nr[1])
    t, wt = roots_legendre(nth)
    th = 0.25 * np.pi * (t + 1.0)
    u, wu = np.cos(th), 0.25 * np.pi * wt * np.sin(th)
    _, ra, rb = distances(R, r, u)
    phi = wf(np.maximum((ra + rb) / R, 1.0), np.clip((ra - rb) / R, -1, 1))
    ls = np.arange(0, lmax + 1, 2)
    P = np.array([eval_legendre(l, u) for l in ls])
    f = (2 * ls[None, :] + 1) * ((phi * wu[None, :]) @ P.T)      # f_l(r) = (2l+1)/2 int_{-1}^1 phi P_l
    # int d3rho phi^2 = sum_l 4 pi/(2l+1) int r^2 f_l^2 dr
    wl = 4 * np.pi / (2 * ls + 1) * np.sum((wr * r * r)[:, None] * f * f, axis=0)
    return ls, wl


# ----------------------------------------------------------------------------- grid choice for the exact states
LEVELS = {  # name: (nR panels, per), (nrho panels, per), ntheta, lmax
    "base": ((16, 10), (24, 12), 128, 64),
    "fine": ((24, 12), (36, 14), 192, 96),
}


def descriptors(M):
    import csv
    fn = os.path.join(ROOT, "paper_data", "mass_scans", "descriptors.csv")
    with open(fn) as fh:
        for row in csv.DictReader(fh):
            if abs(float(row["M"]) - M) <= 1e-6 * M:
                return float(row["R"]), float(row["R2"]), float(row["rXe"])
    raise KeyError(M)


def grid_for(M, level="base", rfac=14.0, sfac=(10.0, 14.0)):
    Rm, R2, rXe = descriptors(M)
    sd = math.sqrt(max(R2 - Rm * Rm, 0.0))
    Rlo = max(0.0, Rm - sfac[0] * sd)
    Rhi = Rm + sfac[1] * sd
    nR, nr, nth, lmax = LEVELS[level]
    return Grid(Rlo, Rhi, rfac * rXe, nR, nr, nth, True, lmax)
