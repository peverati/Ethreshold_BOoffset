"""Product Laguerre basis in perimetric coordinates and matrix construction.

Primitive functions:  F_a(s) = f_{a0}(s0; alpha0) f_{a1}(s1; alpha1) f_{a2}(s2; alpha2)
with index set {a : w0 a0 + w1 a1 + w2 a2 <= N} (default w = (1,1,1), i.e. total
degree <= N; this spans exactly the Hylleraas space
r01^i r02^j r12^k exp(-c01 r01 - c02 r02 - c12 r12), i+j+k <= N, with
alpha0 = c01 + c02, alpha1 = c01 + c12, alpha2 = c02 + c12).

Exchange symmetry (particles 0 <-> 1 is s0 <-> s1, requires alpha0 = alpha1):
    phi_I = n_I (F_a +/- F_a'),  a' = (a1, a0, a2),  a0 >= a1 (a0 > a1 for '-'),
with n_I = 1/2 if a = a' else 1/sqrt(2).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from . import laguerre
from .operators import symmetrize01


@dataclass
class PerimetricBasis:
    N: int
    alphas: Sequence[float]
    symmetry: Optional[str] = None  # None, 'S' (symmetric) or 'A' (antisymmetric) under 0<->1
    weights: Sequence[int] = (1, 1, 1)

    def __post_init__(self):
        self.alphas = tuple(float(a) for a in self.alphas)
        self.weights = tuple(int(w) for w in self.weights)
        if self.symmetry is not None:
            if self.symmetry not in ("S", "A"):
                raise ValueError("symmetry must be None, 'S' or 'A'")
            if self.alphas[0] != self.alphas[1] or self.weights[0] != self.weights[1]:
                raise ValueError("exchange symmetry needs alpha0 == alpha1 and w0 == w1")
        self.full_index = self._full_index()
        self.nmax = tuple(int(self.full_index[:, m].max()) for m in range(3))
        self._build_reduced()

    # ------------------------------------------------------------------ indices
    def _full_index(self):
        w = self.weights
        idx = [
            (a0, a1, a2)
            for a0 in range(self.N // w[0] + 1)
            for a1 in range(self.N // w[1] + 1)
            for a2 in range(self.N // w[2] + 1)
            if w[0] * a0 + w[1] * a1 + w[2] * a2 <= self.N
        ]
        return np.array(idx, dtype=int)

    def _build_reduced(self):
        A = self.full_index
        if self.symmetry is None:
            self.index = A
            self.swapped = None
            self.norms = np.ones(len(A))
            self.sign = 0.0
            return
        if self.symmetry == "S":
            keep = A[:, 0] >= A[:, 1]
            self.sign = 1.0
        else:
            keep = A[:, 0] > A[:, 1]
            self.sign = -1.0
        self.index = A[keep]
        self.swapped = self.index[:, [1, 0, 2]]
        self.norms = np.where(self.index[:, 0] == self.index[:, 1], 0.5, 1.0 / np.sqrt(2.0))

    @property
    def size(self):
        return len(self.index)

    # ------------------------------------------------------------------ 1D data
    def _one_d(self, cache, m, p, da, db):
        key = (m, p, da, db)
        if key not in cache:
            cache[key] = laguerre.one_d_matrix(self.nmax[m], p, da, db, self.alphas[m])
        return cache[key]

    def _raw_matrix(self, terms, A, B, cache):
        """sum_t coef_t prod_m M_m[A_m, B_m]  for row indices A and column indices B."""
        groups = {}
        for coef, p, k, l in terms:
            f = [(p[m], int(k == m), int(l == m)) for m in range(3)]
            groups.setdefault(f[0], {}).setdefault(f[1], []).append((coef, f[2]))
        out = np.zeros((len(A), len(B)))
        r0, c0 = A[:, 0][:, None], B[:, 0][None, :]
        r1, c1 = A[:, 1][:, None], B[:, 1][None, :]
        r2, c2 = A[:, 2][:, None], B[:, 2][None, :]
        for f0, sub in groups.items():
            X0 = self._one_d(cache, 0, *f0)[r0, c0]
            for f1, lst in sub.items():
                C2 = sum(coef * self._one_d(cache, 2, *f2) for coef, f2 in lst)
                X = self._one_d(cache, 1, *f1)[r1, c1] * C2[r2, c2]
                X *= X0
                out += X
        return out

    def matrix(self, terms, cache=None):
        """Matrix of an operator (term list from operators.py) in this basis."""
        if cache is None:
            cache = {}
        if self.symmetry is None:
            Mx = self._raw_matrix(terms, self.index, self.index, cache)
        else:
            t = symmetrize01(terms)
            Mx = self._raw_matrix(t, self.index, self.index, cache)
            Mx += self.sign * self._raw_matrix(t, self.index, self.swapped, cache)
            Mx *= 2.0 * np.outer(self.norms, self.norms)
        return 0.5 * (Mx + Mx.T)

    # --------------------------------------------------------------- coefficients
    def full_tensor(self, c):
        """Coefficients on the primitive product functions as a dense (n0+1, n1+1, n2+1) tensor."""
        T = np.zeros(tuple(n + 1 for n in self.nmax))
        a = self.index
        np.add.at(T, (a[:, 0], a[:, 1], a[:, 2]), self.norms * c)
        if self.symmetry is not None:
            b = self.swapped
            np.add.at(T, (b[:, 0], b[:, 1], b[:, 2]), self.sign * self.norms * c)
        return T

    def describe(self):
        return dict(N=self.N, alphas=list(self.alphas), symmetry=self.symmetry,
                    weights=list(self.weights), size=self.size)
