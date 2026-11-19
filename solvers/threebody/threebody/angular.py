"""Operators for natural-parity states of total angular momentum L (L = 0, 1, 2, ...).

Formulation
-----------
Relative vectors (particle 2 is the "centre"; for {X+, X+, e-} it is the electron,
for helium the nucleus):

    a = x0 - x2   (|a| = r02),      b = x1 - x2   (|b| = r12),      a - b = x0 - x1 (|a-b| = r01).

A state of total angular momentum L, projection M and natural parity (-1)^L can be
written in the vector-coupled (bipolar solid harmonic) form used for two-electron atoms
(A. K. Bhatia and A. Temkin, Rev. Mod. Phys. 36, 1050 (1964); G. W. F. Drake,
Phys. Rev. A 18, 820 (1978); bibliographic data checked online 2026-11-13) as

    Psi_LM = sum_{l1 + l2 = L}  Y^{LM}_{l1 l2}(a, b) f_{l1 l2}(r01, r02, r12),

    Y^{LM}_{l1 l2}(a, b) = |a|^{l1} |b|^{l2} sum_{m1 m2} <l1 m1 l2 m2|L M> Y_{l1 m1}(a^) Y_{l2 m2}(b^),

i.e. L + 1 "components" (bipolar solid harmonics, homogeneous polynomials of degree
l1 in a and l2 in b) times radial functions of the three distances. Here each
component is written as an equivalent *Cartesian* irreducible tensor:

    L = 0 :  A = 1
    L = 1 :  A^(a)_k = a_k ,  A^(b)_k = b_k
    L = 2 :  A^(cd)_kl = (c_k d_l + c_l d_k)/2 - delta_kl (c.d)/3 ,  (cd) in {aa, ab, bb}

(the rank-2 symmetric traceless tensor is the Cartesian form of r^{l1} r^{l2}
[Y_l1 x Y_l2]_2 with l1 + l2 = 2). The map between the 2L+1 spherical components and
the Cartesian components of an irreducible tensor is the same linear map for every
component, so it only contributes a common factor.

Matrix elements. All operators O used here are rotational scalars, so by the
Wigner-Eckart theorem <Psi_LM|O|Psi'_LM'> = delta_MM' <O>_L independent of M, and

    <O>_L  proportional to  sum_K int (A_c,K F)^* O (A_d,K G) dtau

where K runs over all Cartesian components of the tensor. Summed over K the integrand
is invariant under all rotations (and reflections), hence a function of the three
distances only, and the remaining integral over the three Euler angles gives the same
8 pi^2 as for S states (dropped). So the S-state engine (perimetric coordinates,
exact 1D Laguerre integrals) applies unchanged, with new polynomial weights.

Kinetic energy (symmetric, integrated-by-parts form; boundary terms vanish):

    <A_c F | T_i | A_d G> = 1/(2 m_i) sum_K int  grad_i(A_c,K F) . grad_i(A_d,K G)  dtau ,

    grad_i(A F) = F grad_i A + A sum_{j != i} F_{,ij} (x_i - x_j)/r_ij ,

    grad_0 = d/da, grad_1 = d/db, grad_2 = -(d/da + d/db) on functions of (a, b).

The contraction over K and over the Cartesian gradient index is reduced to the three
distances by evaluating the (O(3)-invariant) polynomial in the canonical frame
a = (r02, 0, 0), b = (y, w, 0), y = (r02^2 + r12^2 - r01^2)/(2 r02), w^2 = r12^2 - y^2
(only even powers of w survive; this is asserted). After multiplication by the measure
r01 r02 r12 / 4 every coefficient is a polynomial in the distances (asserted), hence
in the perimetric coordinates: the term format is exactly that of operators.py,
(coef, (p0, p1, p2), k, l) with k, l the perimetric derivative acting on bra / ket.

For L = 0 this reproduces operators.kinetic_terms / overlap_terms / power_terms
exactly (tested).
"""
from __future__ import annotations

import os
import pickle
from functools import lru_cache, wraps

import sympy as sp

from .operators import PAIRS, _DA, _DB, _d_r, _s, _to_terms, combine, third

_a = sp.symbols("a1:4", real=True)
_b = sp.symbols("b1:4", real=True)
_R = {(0, 1): sp.Symbol("R01", positive=True), (0, 2): sp.Symbol("R02", positive=True),
      (1, 2): sp.Symbol("R12", positive=True)}
_y, _w = sp.symbols("y w", real=True)


def components(L: int):
    """Component labels: tuples of vector names, l1 'a's followed by l2 'b's (l1 + l2 = L)."""
    if L < 0:
        raise ValueError("L must be >= 0")
    return tuple(tuple("a" * (L - l2) + "b" * l2) for l2 in range(L + 1))


def swap_component(c):
    """Image of a component under the exchange of particles 0 and 1 (a <-> b).

    For natural parity Y^{LM}_{l1 l2}(b, a) = (-1)^{l1+l2-L} Y^{LM}_{l2 l1}(a, b) = Y^{LM}_{l2 l1}(a, b):
    no phase (Cartesian form: the tensors below are symmetric in their vector arguments).
    """
    return tuple(sorted("b" if x == "a" else "a" for x in c))


def _vec(name):
    return _a if name == "a" else _b


def _dot(u, v):
    return sum(ui * vi for ui, vi in zip(u, v))


def cartesian_tensor(c):
    """Flat list of the Cartesian components of the irreducible tensor of component c."""
    L = len(c)
    if L == 0:
        return [sp.Integer(1)]
    if L == 1:
        return list(_vec(c[0]))
    if L == 2:
        u, v = _vec(c[0]), _vec(c[1])
        uv = _dot(u, v)
        return [sp.Rational(1, 2) * (u[k] * v[l] + u[l] * v[k]) - (sp.Rational(1, 3) * uv if k == l else 0)
                for k in range(3) for l in range(3)]
    raise NotImplementedError("use tensor_components() for L >= 3")


def _harmonic_projection(P, x, n):
    """Harmonic (traceless) part of a homogeneous polynomial P of degree n in the 3-vector x:
    h = sum_j (-1)^j (2n-2j-1)!! / ((2n-1)!! (2j)!!) |x|^(2j) Lap^j P."""
    r2 = sum(xi ** 2 for xi in x)
    lap = lambda f: sum(sp.diff(f, xi, 2) for xi in x)
    out, term = 0, P
    for j in range(n // 2 + 1):
        coef = sp.Integer(-1) ** j * sp.factorial2(2 * n - 2 * j - 1) / (sp.factorial2(2 * n - 1) * sp.factorial2(2 * j))
        out += coef * r2 ** j * term
        term = lap(term)
    return sp.expand(out)


@lru_cache(maxsize=None)
def tensor_components(c):
    """(components, weights) of the rank-L symmetric traceless Cartesian tensor of component c, such that
    sum_K A_c,K A_d,K = sum_k w_k A_c,k A_d,k (and likewise for products of gradients).

    L <= 2: the flat 3^L list of cartesian_tensor() with unit weights (unchanged, so the cached L = 1, 2
    term lists stay valid).  L >= 3: one entry per index multiset (n1, n2, n3), n1 + n2 + n3 = L
    (10 for L = 3, 15 for L = 4), weight = multinomial L!/(n1! n2! n3!).  The tensor is
    T_{i1..iL} = (1/L!) d^L h / dx_{i1} .. dx_{iL}, h = harmonic projection of prod_k (u_k . x),
    i.e. the symmetric traceless part of u_1 (x) ... (x) u_L (Cartesian form of the stretched
    bipolar harmonic [Y_l1(a) Y_l2(b)]_L, l1 + l2 = L)."""
    L = len(c)
    if L <= 2:
        A = cartesian_tensor(c)
        return tuple(A), tuple(sp.Integer(1) for _ in A)
    x = sp.symbols("q1:4", real=True)
    P = sp.Integer(1)
    for name in c:
        P *= _dot(_vec(name), x)
    h = _harmonic_projection(sp.expand(P), x, L)
    comps, weights = [], []
    for n1 in range(L, -1, -1):
        for n2 in range(L - n1, -1, -1):
            n3 = L - n1 - n2
            d = h
            for xi, k in zip(x, (n1, n2, n3)):
                if k:
                    d = sp.diff(d, xi, k)
            comps.append(sp.expand(d / sp.factorial(L)))
            weights.append(sp.factorial(L) / (sp.factorial(n1) * sp.factorial(n2) * sp.factorial(n3)))
    return tuple(comps), tuple(weights)


def _grad_vector(expr, i):
    """Cartesian gradient w.r.t. particle i of a function of (a, b)."""
    da = [sp.diff(expr, x) for x in _a]
    db = [sp.diff(expr, x) for x in _b]
    if i == 0:
        return da
    if i == 1:
        return db
    return [-(p + q) for p, q in zip(da, db)]


def _unit_grad(i, j):
    """grad_i r_ij = (x_i - x_j)/r_ij in terms of a, b."""
    p = (min(i, j), max(i, j))
    if p == (0, 1):
        v = [x - y for x, y in zip(_a, _b)]
    elif p == (0, 2):
        v = list(_a)
    else:
        v = list(_b)
    sign = 1 if i == p[0] else -1
    return [sign * x / _R[p] for x in v]


def _canonical(expr):
    """Reduce an O(3)-invariant polynomial in (a, b) (and 1/R) to a function of the distances."""
    sub = {_a[0]: _R[(0, 2)], _a[1]: 0, _a[2]: 0, _b[0]: _y, _b[1]: _w, _b[2]: 0}
    e = sp.expand(expr.subs(sub))
    P = sp.Poly(e, _w)
    out = 0
    for (k,), coef in P.terms():
        if k % 2:
            if sp.simplify(coef) != 0:
                raise AssertionError("odd power of w survived: expression is not rotation invariant")
            continue
        out += coef * (_R[(1, 2)] ** 2 - _y ** 2) ** (k // 2)
    yv = (_R[(0, 2)] ** 2 + _R[(1, 2)] ** 2 - _R[(0, 1)] ** 2) / (2 * _R[(0, 2)])
    return out.subs(_y, yv)


def _to_perimetric_terms(expr_R):
    """expr in R_ij (must be polynomial after cancel) -> term tuple."""
    num, den = sp.fraction(sp.cancel(sp.together(expr_R)))
    if den.free_symbols & set(_R.values()):
        raise AssertionError(f"non-polynomial integrand, denominator {den}")
    e = (num / den).subs({_R[p]: (_s[p[0]] + _s[p[1]]) / 2 for p in PAIRS})
    return _to_terms(e)


_CACHE_FILE = os.path.join(os.path.dirname(__file__), "_cache", "angular_terms_v1.pkl")
_disk = None


def _disk_cached(fn):
    """Persist the (slow, sympy-derived) term lists in threebody/_cache/ (safe to delete;
    rebuilt in ~5 s for L = 1 and ~40 s for L = 2). Set THREEBODY_NO_TERM_CACHE=1 to bypass."""

    @wraps(fn)
    def wrapper(*args):
        global _disk
        if os.environ.get("THREEBODY_NO_TERM_CACHE"):
            return fn(*args)
        if _disk is None:
            try:
                with open(_CACHE_FILE, "rb") as fh:
                    _disk = pickle.load(fh)
            except Exception:
                _disk = {}
        key = (fn.__name__,) + args
        if key not in _disk:
            _disk[key] = fn(*args)
            try:
                os.makedirs(os.path.dirname(_CACHE_FILE), exist_ok=True)
                old = {}
                try:
                    with open(_CACHE_FILE, "rb") as fh:
                        old = pickle.load(fh)
                except Exception:
                    pass
                old.update(_disk)
                tmp = f"{_CACHE_FILE}.{os.getpid()}.tmp"
                with open(tmp, "wb") as fh:
                    pickle.dump(old, fh)
                os.replace(tmp, _CACHE_FILE)
                _disk = old
            except OSError:  # pragma: no cover - read-only installation
                pass
        return _disk[key]

    return wrapper


def _measure():
    return _R[(0, 1)] * _R[(0, 2)] * _R[(1, 2)] / 4


@lru_cache(maxsize=None)
@_disk_cached
def overlap_terms_L(c, d):
    """Terms of sum_K int (A_c,K F)(A_d,K G) dtau."""
    (Ac, wt), (Ad, _) = tensor_components(c), tensor_components(d)
    expr = sum(w * x * y for w, x, y in zip(wt, Ac, Ad)) * _measure()
    return _to_perimetric_terms(_canonical(expr))


@lru_cache(maxsize=None)
@_disk_cached
def power_terms_L(c, d, i, j, n):
    """Terms of sum_K int (A_c,K F) r_ij^n (A_d,K G) dtau, integer n >= -1."""
    if n < -1:
        raise ValueError("n >= -1 required")
    p = (min(i, j), max(i, j))
    (Ac, wt), (Ad, _) = tensor_components(c), tensor_components(d)
    expr = sum(w * x * y for w, x, y in zip(wt, Ac, Ad)) * _measure() * _R[p] ** n
    return _to_perimetric_terms(_canonical(expr))


@lru_cache(maxsize=None)
@_disk_cached
def kinetic_terms_L(c, d, i):
    """Terms of (1/2) sum_K int grad_i(A_c,K F) . grad_i(A_d,K G) dtau  (multiply by 1/m_i)."""
    (Ac, wt), (Ad, _) = tensor_components(c), tensor_components(d)
    partners = [j for j in (0, 1, 2) if j != i]
    units = {j: _unit_grad(i, j) for j in partners}
    # perimetric marker combinations for d/dr_ij acting on bra (DA) / ket (DB)
    dA = {j: _d_r(min(i, j), max(i, j), _DA) for j in partners}
    dB = {j: _d_r(min(i, j), max(i, j), _DB) for j in partners}
    total = 0
    for wK, AcK, AdK in zip(wt, Ac, Ad):
        gAc, gAd = _grad_vector(AcK, i), _grad_vector(AdK, i)
        for k in range(3):
            bra = gAc[k] + AcK * sum(dA[j] * units[j][k] for j in partners)
            ket = gAd[k] + AdK * sum(dB[j] * units[j][k] for j in partners)
            total += wK * bra * ket
    expr = sp.Rational(1, 2) * sp.expand(total) * _measure()
    # reduce coefficient of each marker monomial separately (keeps sympy fast)
    gens = list(_DA) + list(_DB)
    P = sp.Poly(expr, *gens)
    red = 0
    for monom, coef in P.terms():
        mk = sp.Integer(1)
        for g, e in zip(gens, monom):
            mk *= g ** e
        red += mk * sp.cancel(sp.together(_canonical(coef)))
    return _to_perimetric_terms(red)


def swap_terms(terms):
    """Relabel s0 <-> s1 in a term list (conjugation by the exchange of particles 0 and 1)."""
    perm = {0: 1, 1: 0, 2: 2, -1: -1}
    return tuple((coef, (p[1], p[0], p[2]), perm[k], perm[l]) for coef, p, k, l in terms)


def operator_terms(name, c, d):
    """Term list of a named operator between components c (bra) and d (ket).

    Names: 'S', 'T0', 'T1', 'T2' (per unit 1/m_i), 'V01', 'V02', 'V12' (1/r_ij),
    'T01' = T0 + T1, 'V02+V12', and 'r01', 'r01^2', 'r02', 'r12', 'rXe' = (r02 + r12)/2.
    """
    if name == "S":
        return overlap_terms_L(c, d)
    if name in ("T0", "T1", "T2"):
        return kinetic_terms_L(c, d, int(name[1]))
    if name == "T01":
        return combine((1.0, kinetic_terms_L(c, d, 0)), (1.0, kinetic_terms_L(c, d, 1)))
    if name in ("V01", "V02", "V12"):
        return power_terms_L(c, d, int(name[1]), int(name[2]), -1)
    if name == "V02+V12":
        return combine((1.0, power_terms_L(c, d, 0, 2, -1)), (1.0, power_terms_L(c, d, 1, 2, -1)))
    if name == "r01":
        return power_terms_L(c, d, 0, 1, 1)
    if name == "r01^2":
        return power_terms_L(c, d, 0, 1, 2)
    if name == "1/r01":
        return power_terms_L(c, d, 0, 1, -1)
    if name == "r02":
        return power_terms_L(c, d, 0, 2, 1)
    if name == "r12":
        return power_terms_L(c, d, 1, 2, 1)
    if name == "rXe":
        return combine((0.5, power_terms_L(c, d, 0, 2, 1)), (0.5, power_terms_L(c, d, 1, 2, 1)))
    raise KeyError(name)


__all__ = ["components", "swap_component", "cartesian_tensor", "tensor_components", "overlap_terms_L", "power_terms_L",
           "kinetic_terms_L", "swap_terms", "operator_terms", "third"]
