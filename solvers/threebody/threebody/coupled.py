"""Natural-parity states of total angular momentum L: coupled-component basis and solver.

Wave function (see angular.py for the formulation and conventions):

    Psi_LM = sum_c A_c(a, b) f_c(s0, s1, s2),   c = (l1, l2), l1 + l2 = L,

each radial function f_c expanded in perimetric Laguerre products (one or more
*blocks*, each with its own three scales and index set, as in threshold.py).

Exchange symmetry of identical particles 0 and 1 (P: a <-> b, s0 <-> s1, component
c -> Pc). A symmetry-adapted basis function is

    phi_(c, B, F) = n (A_c F + eps A_Pc PF),   eps = +1 ('S') or -1 ('A'),

with F a primitive of block B and PF its mirror (scales (al1, al0, al2), index
(a1, a0, a2)). One representative c is kept per pair {c, Pc}. For a P-invariant
operator O

    <phi_(c,B,F)| O |phi_(d,C,G)> = 2 n_F n_G ( <A_c F|O|A_d G> + eps <A_c F|O|A_Pd PG> ).

If c = Pc and B is its own mirror (al0 = al1), F and PF give the same function: only
a0 >= a1 is kept (a0 > a1 for eps = -1), n = 1/2 for a0 = a1, else 1/sqrt(2).

Mass-independent matrices: H = sum_i (1/m_i) T_i + sum_ij Z_ij V_ij. For an exchange-
symmetric basis only P-invariant combinations are formed: T01 = T0 + T1, T2, V01,
V02 + V12 (which requires m0 = m1 and Z02 = Z12, guaranteed by `identical`).
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.linalg import eigh

from . import angular
from .threshold import raw_matrix, total_degree_index
from .system import ThreeBodySystem


class IdxBlock:
    """Index set of primitive products with three scales (duck-type compatible with threshold.Block,
    but never truncated automatically)."""

    def __init__(self, alphas, index, name=""):
        self.alphas = tuple(float(a) for a in alphas)
        self.index = np.asarray(index, dtype=int)
        self.name = name

    @property
    def self_mirror(self):
        return self.alphas[0] == self.alphas[1]

    @property
    def nmax(self):
        return tuple(int(self.index[:, m].max()) for m in range(3))

    @property
    def size(self):
        return len(self.index)

    def mirror(self):
        return IdxBlock((self.alphas[1], self.alphas[0], self.alphas[2]), self.index[:, [1, 0, 2]],
                        self.name + "'")

    def subset(self, mask):
        return IdxBlock(self.alphas, self.index[mask], self.name)


def block(N, alphas, weights=(1, 1, 1), name=None):
    """Total-degree block: w0 a0 + w1 a1 + w2 a2 <= N."""
    return IdxBlock(alphas, total_degree_index(N, weights), name or f"N{N}")


SYM_OPS = ("S", "T01", "T2", "V01", "V02+V12")
GEN_OPS = ("S", "T0", "T1", "T2", "V01", "V02", "V12")


class CoupledBasis:
    """Basis for natural-parity total angular momentum L.

    L        : total angular momentum (0, 1, 2)
    blocks   : list of IdxBlock (or a single one); used for every component
    symmetry : 'S', 'A' (exchange of particles 0 and 1) or None (no symmetry)
    """

    def __init__(self, L: int, blocks, symmetry: Optional[str] = "S", ops=None, comp_blocks=None):
        self.L = int(L)
        self.blocks = [blocks] if isinstance(blocks, IdxBlock) else list(blocks)
        if symmetry not in ("S", "A", None):
            raise ValueError("symmetry must be 'S', 'A' or None")
        self.symmetry = symmetry
        self.eps = {"S": 1.0, "A": -1.0, None: 0.0}[symmetry]
        comps = angular.components(self.L)
        if symmetry is None:
            self.channels = list(comps)
        else:
            self.channels = [c for c in comps if c <= angular.swap_component(c)]
        # optional per-channel block lists (default: all blocks for every channel)
        self.comp_blocks = comp_blocks or {c: self.blocks for c in self.channels}
        self.units = []  # list of (component, IdxBlock (possibly truncated), norms)
        for c in self.channels:
            for B in self.comp_blocks[c]:
                self.units.append(self._make_unit(c, B))
        self.ops = tuple(ops) if ops is not None else (SYM_OPS if symmetry else GEN_OPS)
        self.offsets = np.cumsum([0] + [u[1].size for u in self.units])
        t0 = time.perf_counter()
        self.mats = {o: self._assemble(o) for o in self.ops}
        self.build_time = time.perf_counter() - t0
        self._proj = None

    # ------------------------------------------------------------------ units
    def _make_unit(self, c, B: IdxBlock):
        if self.symmetry is not None and angular.swap_component(c) == c and B.self_mirror:
            a = B.index
            keep = a[:, 0] >= a[:, 1] if self.eps > 0 else a[:, 0] > a[:, 1]
            Bt = B.subset(keep)
            norms = np.where(Bt.index[:, 0] == Bt.index[:, 1], 0.5, 1.0 / np.sqrt(2.0))
            return (c, Bt, norms)
        n = 1.0 / np.sqrt(2.0) if self.symmetry is not None else 1.0
        return (c, B, np.full(B.size, n))

    @property
    def size(self):
        return int(self.offsets[-1])

    def describe(self):
        return dict(L=self.L, symmetry=self.symmetry, size=self.size,
                    units=[dict(component="".join(c) or "s", alphas=list(B.alphas), size=B.size, nmax=list(B.nmax),
                                name=B.name) for c, B, _ in self.units])

    # --------------------------------------------------------------- assembly
    def _assemble(self, op):
        n = self.size
        out = np.zeros((n, n))
        cache = {}
        off = self.offsets
        for I, (c, B, nb) in enumerate(self.units):
            for J in range(I, len(self.units)):
                d, C, nc = self.units[J]
                blk = raw_matrix(angular.operator_terms(op, c, d), B, C, cache)
                if self.symmetry is not None:
                    Pd = angular.swap_component(d)
                    blk = blk + self.eps * raw_matrix(angular.operator_terms(op, c, Pd), B, C.mirror(), cache)
                    blk *= 2.0 * np.outer(nb, nc)
                out[off[I]:off[I + 1], off[J]:off[J + 1]] = blk
                if J != I:
                    out[off[J]:off[J + 1], off[I]:off[I + 1]] = blk.T
        return 0.5 * (out + out.T)

    # ------------------------------------------------------------ hamiltonian
    def hamiltonian(self, system: ThreeBodySystem):
        im = system.inverse_masses
        z = system.coupling_dict()
        m = self.mats
        if self.symmetry is not None:
            if system.identical is None:
                raise ValueError("symmetry-adapted basis needs identical particles 0 and 1")
            H = im[0] * m["T01"]
            if im[2] != 0.0:
                H = H + im[2] * m["T2"]
            H = H + z[(0, 1)] * m["V01"] + z[(0, 2)] * m["V02+V12"]
            return H
        H = sum(im[i] * m[f"T{i}"] for i in range(3) if im[i] != 0.0)
        for (i, j), zz in z.items():
            if zz != 0.0:
                H = H + zz * m[f"V{i}{j}"]
        return H

    def project(self, tol=1e-12):
        """Canonical orthogonalization (S scaled to unit diagonal; drop eigenvalues < tol * max)."""
        if self._proj is not None and self._proj["tol"] == tol:
            return self._proj
        S = self.mats["S"]
        dd = 1.0 / np.sqrt(np.diag(S))
        lam, U = np.linalg.eigh(S * np.outer(dd, dd))
        keep = lam > tol * lam[-1]
        X = (dd[:, None] * U[:, keep]) / np.sqrt(lam[keep])
        self._proj = dict(tol=tol, X=X, kept=int(keep.sum()), dropped=int((~keep).sum()),
                          lam_min=float(lam[0]), lam_max=float(lam[-1]), P={})
        return self._proj

    def projected(self, op, tol=1e-12):
        pr = self.project(tol)
        if op not in pr["P"]:
            pr["P"][op] = pr["X"].T @ self.mats[op] @ pr["X"]
        return pr["P"][op]

    def hamiltonian_projected(self, system, tol=1e-12):
        """Projected H; for exchange-symmetric bases built from the cached projected matrices."""
        pr = self.project(tol)
        if self.symmetry is not None:
            im = system.inverse_masses
            z = system.coupling_dict()
            H = im[0] * self.projected("T01", tol)
            if im[2] != 0.0:
                H = H + im[2] * self.projected("T2", tol)
            return H + z[(0, 1)] * self.projected("V01", tol) + z[(0, 2)] * self.projected("V02+V12", tol)
        return pr["X"].T @ self.hamiltonian(system) @ pr["X"]

    def solve(self, system, nstates=3, tol=1e-12, vectors=True):
        H = self.hamiltonian_projected(system, tol)
        k = min(nstates, H.shape[0])
        if not vectors:
            return eigh(H, eigvals_only=True, subset_by_index=[0, k - 1]), None
        e, U = eigh(H, subset_by_index=[0, k - 1])
        return e, self.project(tol)["X"] @ U  # coefficients in the (non-orthogonal) basis, c^T S c = 1

    def expectation(self, op, C):
        """<psi|op|psi> for normalized coefficient vector(s) C (columns); op must be assembled or assemblable."""
        if op not in self.mats:
            self.mats[op] = self._assemble(op)
        C = np.atleast_2d(C.T).T
        return np.einsum("ia,ij,ja->a", C, self.mats[op], C)


def coupled_threshold(system: ThreeBodySystem):
    """Continuum edge for natural-parity states of any L (see README): lowest two-body
    threshold, the atom in its 1s state and relative angular momentum lambda = L."""
    return system.threshold()


def solve_L(system, L, N, alphas, symmetry="S", weights=(1, 1, 1), nstates=3, tol=1e-12):
    """Convenience: single total-degree block for every component."""
    b = CoupledBasis(L, block(N, alphas, weights), symmetry)
    e, C = b.solve(system, nstates, tol)
    return e, C, b
