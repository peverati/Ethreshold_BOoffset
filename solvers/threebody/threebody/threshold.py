"""Multi-scale perimetric basis for near-threshold S states of {X+, X+, e-}.

Why: close to the X + (Xe) threshold an excited state is an ion weakly bound to a
neutral atom.  Its wave function is compact in the atom coordinates (one X with the
electron) and extremely diffuse in the ion-atom distance R.  The package basis has
one scale per perimetric coordinate and alpha0 = alpha1 (needed for exchange
symmetry), so it cannot be compact in s0 and diffuse in s1 at the same time.

Here the basis is a union of *blocks*.  A block is a set of product functions

    F_a(s) = f_a0(s0; al0) f_a1(s1; al1) f_a2(s2; al2),  a in index set,

with its own three scales (al0 != al1 allowed).  Exchange symmetry (0 <-> 1, i.e.
s0 <-> s1, P) is imposed by using phi = F + P F, where P F has scales
(al1, al0, al2) and index (a1, a0, a2).  For an exchange-symmetric operator O

    <F + PF | O | G + PG> = 2 (<F|O|G> + <F|O|PG>),

so only "direct" and "mirror" raw blocks are needed.  A self-mirror block
(al0 = al1) keeps a0 >= a1 only (otherwise phi would be duplicated).

Typical union:  compact block (the usual total-degree basis, scales optimized for
the excited root) + diffuse block with scales (mu_a, al_d, mu_a), small degree in
(s0, s2) (atom + its polarization) and large degree in s1 (s1 ~ 2R when X1 is far
from the atom (X0 e), so exp(-al_d s1/2) ~ exp(-al_d R)).

All 1D integrals are exact: same scale -> laguerre.one_d_matrix; different scales
-> Gauss-Laguerre (integrand is polynomial x exponential, see
laguerre.mixed_moment_matrix), derivatives via d f_n/ds = al sum_k D_kn f_k.

H(M) = (1/M) T_X + T_e + V_01 - V_Xe with mass-independent matrices, so a mass
scan at fixed basis only re-diagonalizes a small projected problem.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import eigh

from . import laguerre
from .operators import combine, kinetic_terms, overlap_terms, power_terms
from .system import family_threshold


# ---------------------------------------------------------------- index sets
def total_degree_index(N, weights=(1, 1, 1)):
    w = weights
    return np.array([(a0, a1, a2)
                     for a0 in range(N // w[0] + 1) for a1 in range(N // w[1] + 1) for a2 in range(N // w[2] + 1)
                     if w[0] * a0 + w[1] * a1 + w[2] * a2 <= N], dtype=int)


def diffuse_index(n_c, n_d, n_cross=None):
    """a0 + a2 <= n_c (atom coordinates), a1 <= n_d (ion-atom coordinate s1).

    n_cross (optional) limits the coupling: a0 + a2 <= n_cross when a1 > n_d//2
    is NOT used by default; kept simple (full product)."""
    return np.array([(a0, a1, a2) for a0 in range(n_c + 1) for a2 in range(n_c + 1 - a0)
                     for a1 in range(n_d + 1)], dtype=int)


@dataclass
class Block:
    alphas: tuple
    index: np.ndarray
    name: str = ""

    def __post_init__(self):
        self.alphas = tuple(float(a) for a in self.alphas)
        self.index = np.asarray(self.index, dtype=int)
        if self.self_mirror:
            self.index = self.index[self.index[:, 0] >= self.index[:, 1]]

    @property
    def self_mirror(self):
        return self.alphas[0] == self.alphas[1]

    @property
    def nmax(self):
        return tuple(int(self.index[:, m].max()) for m in range(3))

    def mirror(self):
        b = Block.__new__(Block)
        b.alphas = (self.alphas[1], self.alphas[0], self.alphas[2])
        b.index = self.index[:, [1, 0, 2]]
        b.name = self.name + "'"
        return b

    @property
    def size(self):
        return len(self.index)


def compact_block(N, alpha_x, alpha_e, weights=(1, 1, 1)):
    """Usual symmetric total-degree block (same function space as PerimetricBasis(N, 'S'))."""
    return Block((alpha_x, alpha_x, alpha_e), total_degree_index(N, weights), name=f"C{N}")


def diffuse_block(alpha_atom, alpha_d, n_c, n_d):
    return Block((alpha_atom, alpha_d, alpha_atom), diffuse_index(n_c, n_d), name=f"D{n_c}x{n_d}@{alpha_d:g}")


# ---------------------------------------------------------------- 1D integrals
def one_d_mixed(alpha_a, na, alpha_b, nb, p, da, db):
    """int (d^da f_a(s;alpha_a)) (d^db f_b(s;alpha_b)) s^p ds, a <= na, b <= nb (exact)."""
    if alpha_a == alpha_b:
        n = max(na, nb)
        return laguerre.one_d_matrix(n, p, da, db, alpha_a)[: na + 1, : nb + 1]
    S = laguerre.mixed_moment_matrix(na, alpha_a, nb, alpha_b, p)
    if da:
        S = alpha_a * (laguerre.derivative_matrix(na + 1).T @ S)
    if db:
        S = alpha_b * (S @ laguerre.derivative_matrix(nb + 1))
    return S


def raw_matrix(terms, A: Block, B: Block, cache):
    """sum_t coef_t prod_m int s_m^p (d f_{A,m}) (d f_{B,m}) ds_m between two (unsymmetrized) blocks."""
    na, nb = A.nmax, B.nmax

    def od(m, p, da, db):
        key = (A.alphas[m], na[m], B.alphas[m], nb[m], p, da, db)
        if key not in cache:
            cache[key] = one_d_mixed(A.alphas[m], na[m], B.alphas[m], nb[m], p, da, db)
        return cache[key]

    groups = {}
    for coef, p, k, l in terms:
        f = [(p[m], int(k == m), int(l == m)) for m in range(3)]
        groups.setdefault(f[0], {}).setdefault(f[1], []).append((coef, f[2]))
    out = np.zeros((A.size, B.size))
    r0, c0 = A.index[:, 0][:, None], B.index[:, 0][None, :]
    r1, c1 = A.index[:, 1][:, None], B.index[:, 1][None, :]
    r2, c2 = A.index[:, 2][:, None], B.index[:, 2][None, :]
    for f0, sub in groups.items():
        X0 = od(0, *f0)[r0, c0]
        for f1, lst in sub.items():
            C2 = sum(coef * od(2, *f2) for coef, f2 in lst)
            X = od(1, *f1)[r1, c1] * C2[r2, c2]
            X *= X0
            out += X
    return out


# ---------------------------------------------------------------- operators of the family
def family_operators():
    """Exchange-symmetric, mass-independent operator term lists."""
    k = [kinetic_terms(i) for i in range(3)]
    return {
        "S": overlap_terms(),
        "TX": combine((1.0, k[0]), (1.0, k[1])),          # times 1/M
        "Te": k[2],                                       # electron, mass 1
        "V01": power_terms(0, 1, -1),                     # +1/r01
        "VXe": combine((1.0, power_terms(0, 2, -1)), (1.0, power_terms(1, 2, -1))),  # times -1
        "r01": power_terms(0, 1, 1),
        "r01^2": power_terms(0, 1, 2),
        "rXe": combine((0.5, power_terms(0, 2, 1)), (0.5, power_terms(1, 2, 1))),
    }


class MultiBlockBasis:
    """Exchange-symmetric union of blocks; mass-independent matrices of the {M,M,1} family."""

    def __init__(self, blocks, ops=("S", "TX", "Te", "V01", "VXe", "r01", "r01^2", "rXe"), sign=+1.0):
        self.blocks = list(blocks)
        self.sign = sign  # +1: symmetric (phi = F + PF); -1: antisymmetric
        t0 = time.perf_counter()
        allops = family_operators()
        self.mats = {o: self._assemble(allops[o]) for o in ops}
        self.build_time = time.perf_counter() - t0
        self._proj = None

    @property
    def size(self):
        return sum(b.size for b in self.blocks)

    def describe(self):
        return [dict(name=b.name, alphas=list(b.alphas), size=b.size, nmax=list(b.nmax)) for b in self.blocks]

    def _assemble(self, terms):
        n = self.size
        out = np.zeros((n, n))
        off = np.cumsum([0] + [b.size for b in self.blocks])
        cache = {}
        for I, A in enumerate(self.blocks):
            for J in range(I, len(self.blocks)):
                B = self.blocks[J]
                blk = raw_matrix(terms, A, B, cache) + self.sign * raw_matrix(terms, A, B.mirror(), cache)
                out[off[I]:off[I + 1], off[J]:off[J + 1]] = 2.0 * blk
                if J != I:
                    out[off[J]:off[J + 1], off[I]:off[I + 1]] = 2.0 * blk.T
        return 0.5 * (out + out.T)

    # ------------------------------------------------------------ projection
    def project(self, tol=1e-12):
        """Canonical orthogonalization (S scaled to unit diagonal, drop eig < tol*max)."""
        if self._proj is not None and self._proj["tol"] == tol:
            return self._proj
        S = self.mats["S"]
        d = 1.0 / np.sqrt(np.diag(S))
        lam, U = np.linalg.eigh(S * np.outer(d, d))
        keep = lam > tol * lam[-1]
        X = (d[:, None] * U[:, keep]) / np.sqrt(lam[keep])
        P = {o: X.T @ m @ X for o, m in self.mats.items() if o != "S"}
        self._proj = dict(tol=tol, X=X, P=P, kept=int(keep.sum()), dropped=int((~keep).sum()),
                          lam_min=float(lam[0]), lam_max=float(lam[-1]))
        return self._proj

    def hamiltonian_projected(self, M, tol=1e-12):
        P = self.project(tol)["P"]
        return P["TX"] / M + P["Te"] + P["V01"] - P["VXe"]

    def energies(self, M, nstates=3, tol=1e-12):
        H = self.hamiltonian_projected(M, tol)
        return eigh(H, eigvals_only=True, subset_by_index=[0, nstates - 1])

    def solve(self, M, nstates=3, tol=1e-12):
        H = self.hamiltonian_projected(M, tol)
        e, U = eigh(H, subset_by_index=[0, nstates - 1])
        pr = self.project(tol)
        return e, U, pr

    def observables(self, M, state=1, tol=1e-12):
        """Expectation values (normalized) of r01, r01^2, r_Xe and dE/dM (Hellmann-Feynman, fixed basis)."""
        e, U, pr = self.solve(M, state + 1, tol)
        u = U[:, state]
        P = pr["P"]
        val = lambda o: float(u @ P[o] @ u)
        TX = val("TX")
        return dict(E=float(e[state]), r01=val("r01"), r01_2=val("r01^2"), rXe=val("rXe"),
                    T_X=TX / M, dE_dM=-TX / M ** 2)

    def critical_mass(self, M_lo, M_hi, state=1, tol=1e-12, xtol=1e-10):
        """Mass at which root `state` of this fixed basis crosses E_th(M) (Brent).

        Since root `state` is a rigorous upper bound to the exact level, this M* is
        a rigorous *upper bound* to the exact critical mass (up to rounding)."""
        from scipy.optimize import brentq

        f = lambda M: self.energies(M, state + 1, tol)[state] - family_threshold(M)
        flo, fhi = f(M_lo), f(M_hi)
        if not (flo > 0 > fhi):
            raise ValueError(f"no sign change: f({M_lo})={flo:.3e}, f({M_hi})={fhi:.3e}")
        return brentq(f, M_lo, M_hi, xtol=xtol, rtol=1e-14)


def threshold_reduced_masses(M):
    """(atom reduced mass M/(M+1), ion-atom reduced mass M(M+1)/(2M+1))."""
    return M / (M + 1.0), M * (M + 1.0) / (2.0 * M + 1.0)
