"""Operators of the three-body problem (S states) in perimetric coordinates.

Conventions
-----------
Particles are labelled 0, 1, 2; r_ij is the distance between i and j.
Perimetric coordinates (one per vertex of the particle triangle):

    s_k = r_kj + r_kl - r_jl      ({k, j, l} = {0, 1, 2}),   s_k in [0, inf)

so that r_ij = (s_i + s_j) / 2 and the S-state volume element is

    dtau = 8 pi^2 r01 r02 r12 dr01 dr02 dr12 = 8 pi^2 (r01 r02 r12 / 4) ds0 ds1 ds2 .

The constant 8 pi^2 is dropped everywhere (it cancels in every ratio).

Every operator matrix element <phi_a | O | phi_b> is reduced to a sum of terms

    coef * int  s0^p0 s1^p1 s2^p2  (d_k phi_a) (d_l phi_b)  ds

where d_k = d/ds_k (k = -1 means "no derivative").  A term is stored as
(coef, (p0, p1, p2), k, l).  Because the basis is a product of 1D functions,
each term is a product of three exact 1D integrals (see laguerre.py).

Kinetic energy (derivation, verified symbolically in tests/test_hamiltonian.py):
for f = f(r01, r02, r12) and particle i with partners j, k

    grad_i f = f_{,ij} rhat_ij + f_{,ik} rhat_ik ,  rhat_ij . rhat_ik = cos_i
             = (r_ij^2 + r_ik^2 - r_jk^2) / (2 r_ij r_ik)

so the (integrated-by-parts, symmetric) kinetic matrix element is

    <a|T|b> = sum_i 1/(2 m_i) int [a_ij b_ij + a_ik b_ik + cos_i (a_ij b_ik + a_ik b_ij)] dtau ,

and cos_i * (r01 r02 r12) = (r_ij^2 + r_ik^2 - r_jk^2) r_jk / 2 is a polynomial:
no singular integrand ever appears.  In perimetric coordinates
d/dr_ij = d/ds_i + d/ds_j - d/ds_k.
"""
from __future__ import annotations

from functools import lru_cache

import sympy as sp

PAIRS = ((0, 1), (0, 2), (1, 2))

_s = sp.symbols("s0 s1 s2", nonnegative=True)
_DA = sp.symbols("DA0 DA1 DA2")  # markers for d/ds_k acting on the bra
_DB = sp.symbols("DB0 DB1 DB2")  # markers for d/ds_k acting on the ket


def third(i: int, j: int) -> int:
    return 3 - i - j


def r(i: int, j: int):
    return (_s[i] + _s[j]) / 2


def _d_r(i: int, j: int, markers):
    k = third(i, j)
    return markers[i] + markers[j] - markers[k]


def _weight():
    """r01 r02 r12 times the Jacobian 1/4 of (r01,r02,r12) -> (s0,s1,s2)."""
    return r(0, 1) * r(0, 2) * r(1, 2) / 4


def _to_terms(expr) -> tuple:
    """Convert a sympy expression (polynomial in s, at most bilinear in DA/DB) to terms."""
    gens = list(_s) + list(_DA) + list(_DB)
    poly = sp.Poly(sp.expand(expr), *gens)
    terms = []
    for monom, coef in poly.terms():
        p = monom[0:3]
        da = monom[3:6]
        db = monom[6:9]
        if sum(da) > 1 or sum(db) > 1:
            raise ValueError("operator is not bilinear in first derivatives")
        k = da.index(1) if sum(da) == 1 else -1
        l = db.index(1) if sum(db) == 1 else -1
        terms.append((float(coef), tuple(int(x) for x in p), k, l))
    return tuple(terms)


@lru_cache(maxsize=None)
def overlap_terms():
    return _to_terms(_weight())


@lru_cache(maxsize=None)
def kinetic_terms(i: int):
    """Terms of (1/2) int grad_i a . grad_i b dtau (multiply by 1/m_i)."""
    j, k = [x for x in (0, 1, 2) if x != i]
    a_ij, a_ik = _d_r(i, j, _DA), _d_r(i, k, _DA)
    b_ij, b_ik = _d_r(i, j, _DB), _d_r(i, k, _DB)
    rij, rik, rjk = r(i, j), r(i, k), r(j, k)
    w = _weight()
    cos_w = (rij**2 + rik**2 - rjk**2) * rjk / 2 / 4  # cos_i * r01 r02 r12 / 4
    expr = sp.Rational(1, 2) * ((a_ij * b_ij + a_ik * b_ik) * w + (a_ij * b_ik + a_ik * b_ij) * cos_w)
    return _to_terms(expr)


@lru_cache(maxsize=None)
def power_terms(i: int, j: int, n: int):
    """Terms of r_ij^n for integer n >= -1 (n = -1 is the Coulomb operator)."""
    if n < -1:
        raise ValueError("r_ij^n with n < -1 is not polynomial; use pair-distribution quadrature")
    i, j = min(i, j), max(i, j)
    k = third(i, j)
    if n == -1:
        expr = r(i, k) * r(j, k) / 4
    else:
        expr = _weight() * r(i, j) ** n
    return _to_terms(expr)


def swap01(terms):
    """Terms of the operator conjugated by the exchange of particles 0 and 1."""
    perm = {0: 1, 1: 0, 2: 2, -1: -1}
    out = []
    for coef, p, k, l in terms:
        out.append((coef, (p[1], p[0], p[2]), perm[k], perm[l]))
    return tuple(out)


def combine(*weighted):
    """Linear combination: combine((c1, terms1), (c2, terms2), ...)."""
    acc = {}
    for c, terms in weighted:
        if c == 0:
            continue
        for coef, p, k, l in terms:
            key = (p, k, l)
            acc[key] = acc.get(key, 0.0) + c * coef
    return tuple((v, key[0], key[1], key[2]) for key, v in acc.items() if v != 0.0)


def symmetrize01(terms):
    return combine((0.5, terms), (0.5, swap01(terms)))
