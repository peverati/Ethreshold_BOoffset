"""Verification of the L > 0 (natural parity) operators and coupled basis.

1. L = 0 term lists from angular.py are identical to operators.py.
2. Pointwise: the term-list kinetic integrand equals a brute-force Cartesian
   evaluation of sum_K grad_i(A_cK F) . grad_i(A_dK G) (sympy, 9 coordinates) at
   random configurations, for every particle and every component pair (L = 1, 2).
3. Matrix elements: engine (first-derivative form) versus Gauss-Laguerre quadrature
   of the second-order form <A_c F| -1/2 Lap_i |A_d G> with a Cartesian Laplacian
   (no symmetry, three different scales), plus overlap and Coulomb.
4. Exchange symmetry: S + A spectra equal the unsymmetrized spectrum.
5. Hermiticity of the assembled matrices.
"""
import numpy as np
import pytest
import sympy as sp
from scipy.special import roots_laguerre

from threebody import angular, helium_like, symmetric_family
from threebody.coupled import CoupledBasis, block
from threebody.operators import combine, kinetic_terms, overlap_terms, power_terms, PAIRS

X = sp.symbols("x0:9", real=True)
POS = [sp.Matrix(X[3 * n:3 * n + 3]) for n in range(3)]
VA = POS[0] - POS[2]
VB = POS[1] - POS[2]
DIST = {p: sp.sqrt(((POS[p[0]] - POS[p[1]]).T * (POS[p[0]] - POS[p[1]]))[0]) for p in PAIRS}
S_OF_X = (DIST[(0, 1)] + DIST[(0, 2)] - DIST[(1, 2)], DIST[(0, 1)] + DIST[(1, 2)] - DIST[(0, 2)],
          DIST[(0, 2)] + DIST[(1, 2)] - DIST[(0, 1)])
SS = sp.symbols("s0:3", positive=True)


def _tensor_cart(c):
    """Cartesian components of the irreducible tensor of component c in terms of x (independent code)."""
    vec = {"a": VA, "b": VB}
    if len(c) == 0:
        return [sp.Integer(1)]
    if len(c) == 1:
        return [vec[c[0]][k] for k in range(3)]
    if len(c) == 2:
        u, v = vec[c[0]], vec[c[1]]
        uv = (u.T * v)[0]
        return [(u[k] * v[l] + u[l] * v[k]) / 2 - (uv / 3 if k == l else 0) for k in range(3) for l in range(3)]
    return _traceless_flat([vec[n] for n in c])


def _sym_flat(vs):
    """Symmetrized outer product of the vectors vs as a dict over index tuples (independent of angular.py)."""
    import itertools
    L = len(vs)
    out = {}
    for idx in itertools.product(range(3), repeat=L):
        out[idx] = sum(sp.prod([vs[p[k]][idx[k]] for k in range(L)]) for p in itertools.permutations(range(L))) \
            / sp.factorial(L)
    return out


def _traceless_flat(vs):
    """Rank 3 / 4 symmetric traceless part by the explicit detracer formulas (independent of angular.py):
    L=3: T_ijk = S_ijk - (1/5)(d_ij S_llk + d_ik S_ljl + d_jk S_ill)
    L=4: T_ijkl = S_ijkl - (1/7)(sum of 6 d S..mm) + (1/35)(d_ij d_kl + d_ik d_jl + d_il d_jk) S_mmnn."""
    import itertools
    L = len(vs)
    S = _sym_flat(vs)
    d = lambda i, j: 1 if i == j else 0
    if L == 3:
        tr = {k: sum(S[(m, m, k)] for m in range(3)) for k in range(3)}
        return [S[(i, j, k)] - sp.Rational(1, 5) * (d(i, j) * tr[k] + d(i, k) * tr[j] + d(j, k) * tr[i])
                for i, j, k in itertools.product(range(3), repeat=3)]
    if L == 4:
        tr = {(k, l): sum(S[(m, m, k, l)] for m in range(3)) for k in range(3) for l in range(3)}
        tt = sum(tr[(n, n)] for n in range(3))
        out = []
        for i, j, k, l in itertools.product(range(3), repeat=4):
            t2 = (d(i, j) * tr[(k, l)] + d(i, k) * tr[(j, l)] + d(i, l) * tr[(j, k)] + d(j, k) * tr[(i, l)]
                  + d(j, l) * tr[(i, k)] + d(k, l) * tr[(i, j)])
            t4 = (d(i, j) * d(k, l) + d(i, k) * d(j, l) + d(i, l) * d(j, k)) * tt
            out.append(S[(i, j, k, l)] - sp.Rational(1, 7) * t2 + sp.Rational(1, 35) * t4)
        return out
    raise NotImplementedError


def _same(t1, t2):
    d = combine((1.0, t1), (-1.0, t2))
    return max([abs(x[0]) for x in d] + [0.0])


def test_L0_reduces_to_S_state_operators():
    for i in range(3):
        assert _same(angular.kinetic_terms_L((), (), i), kinetic_terms(i)) < 1e-14
    assert _same(angular.overlap_terms_L((), ()), overlap_terms()) < 1e-14
    for p in PAIRS:
        assert _same(angular.power_terms_L((), (), p[0], p[1], -1), power_terms(p[0], p[1], -1)) < 1e-14


def _eval_terms(terms, F, G, s):
    """sum coef s^p (d_k F)(d_l G) at the point s (F, G sympy expressions in SS)."""
    sub = dict(zip(SS, s))
    dF = {-1: float(F.subs(sub)), **{k: float(sp.diff(F, SS[k]).subs(sub)) for k in range(3)}}
    dG = {-1: float(G.subs(sub)), **{k: float(sp.diff(G, SS[k]).subs(sub)) for k in range(3)}}
    return sum(coef * s[0] ** p[0] * s[1] ** p[1] * s[2] ** p[2] * dF[k] * dG[l] for coef, p, k, l in terms)


@pytest.mark.parametrize("L", [1, 2, 3])
def test_kinetic_integrand_pointwise_against_cartesian(L):
    F = (1 + sp.Rational(3, 10) * SS[0] + sp.Rational(1, 5) * SS[1] * SS[2] + SS[0] ** 2 / 10) * sp.exp(
        -sp.Rational(7, 10) * SS[0] - sp.Rational(2, 5) * SS[1] - sp.Rational(9, 10) * SS[2])
    G = (1 + SS[2] / 3 + SS[0] * SS[1] / 7) * sp.exp(-SS[0] / 2 - sp.Rational(6, 5) * SS[1] - SS[2] / 4)
    Fx = F.subs(dict(zip(SS, S_OF_X)))
    Gx = G.subs(dict(zip(SS, S_OF_X)))
    rng = np.random.default_rng(7 + L)
    comps = angular.components(L)
    pts = [rng.normal(size=9) for _ in range(2)]
    for c in comps:
        for d in comps:
            Ac, Ad = _tensor_cart(c), _tensor_cart(d)
            for i in range(3):
                terms = angular.kinetic_terms_L(c, d, i)
                xs = X[3 * i:3 * i + 3]
                expr = sum(sum(sp.diff(a * Fx, x) * sp.diff(b * Gx, x) for x in xs) for a, b in zip(Ac, Ad))
                fn = sp.lambdify(X, expr, "numpy")
                for x in pts:
                    sub = dict(zip(X, x))
                    r = {p: float(DIST[p].subs(sub)) for p in PAIRS}
                    s = [r[(0, 1)] + r[(0, 2)] - r[(1, 2)], r[(0, 1)] + r[(1, 2)] - r[(0, 2)],
                         r[(0, 2)] + r[(1, 2)] - r[(0, 1)]]
                    meas = r[(0, 1)] * r[(0, 2)] * r[(1, 2)] / 4
                    lhs = 2.0 * _eval_terms(terms, F, G, s) / meas
                    rhs = float(fn(*x))
                    assert abs(lhs - rhs) < 1e-10 * max(1.0, abs(rhs)), (c, d, i, lhs, rhs)


def _quadrature_check(L, alphas, N, nq=22):
    basis = CoupledBasis(L, block(N, alphas), None, ops=("S", "T0", "T1", "T2", "V01", "V02", "V12"))
    # primitive functions in s
    funcs = []
    for c, B, _ in basis.units:
        for a in B.index:
            f = 1
            for m in range(3):
                f *= sp.assoc_laguerre(int(a[m]), 0, alphas[m] * SS[m]) * sp.exp(-alphas[m] * SS[m] / 2)
            funcs.append((c, f))
    x, w = roots_laguerre(nq)
    grids = [(x / al, w * np.exp(x) / al) for al in alphas]
    S0, S1, S2 = np.meshgrid(grids[0][0], grids[1][0], grids[2][0], indexing="ij")
    W = np.einsum("i,j,k->ijk", grids[0][1], grids[1][1], grids[2][1])
    R01, R02, R12 = (S0 + S1) / 2, (S0 + S2) / 2, (S1 + S2) / 2
    meas = R01 * R02 * R12 / 4
    # canonical configuration: x2 = 0, x0 = a = (r02, 0, 0), x1 = b = (y, w, 0)
    yv = (R02 ** 2 + R12 ** 2 - R01 ** 2) / (2 * R02)
    wv = np.sqrt(np.maximum(R12 ** 2 - yv ** 2, 0.0))
    zero = np.zeros_like(R01)
    coords = (R02, zero, zero, yv, wv, zero, zero, zero, zero)
    n = len(funcs)
    ops = {}
    psi = []  # list of per-K value arrays for each function
    for c, f in funcs:
        A = _tensor_cart(c)
        fx = f.subs(dict(zip(SS, S_OF_X)))
        psi.append([a * fx for a in A])
    vals = [[sp.lambdify(X, comp, "numpy")(*coords) * np.ones_like(R01) for comp in p] for p in psi]
    S_ref = np.array([[np.sum(W * meas * sum(u * v for u, v in zip(vals[p], vals[q]))) for q in range(n)]
                      for p in range(n)])
    ops["S"] = S_ref
    for (i, j) in PAIRS:
        r = {(0, 1): R01, (0, 2): R02, (1, 2): R12}[(i, j)]
        ops[f"V{i}{j}"] = np.array([[np.sum(W * meas / r * sum(u * v for u, v in zip(vals[p], vals[q])))
                                     for q in range(n)] for p in range(n)])
    for i in range(3):
        xs = X[3 * i:3 * i + 3]
        lap = [[sp.lambdify(X, -sum(sp.diff(comp, x, 2) for x in xs) / 2, "numpy")(*coords) * np.ones_like(R01)
                for comp in p] for p in psi]
        ops[f"T{i}"] = np.array([[np.sum(W * meas * sum(u * v for u, v in zip(vals[p], lap[q])))
                                  for q in range(n)] for p in range(n)])
    return basis, ops


@pytest.mark.parametrize("L", [1, 2])
def test_matrix_elements_against_quadrature_second_order_form(L):
    alphas = (0.9, 1.4, 1.7)
    # L = 2: one primitive per component (all 9 component pairs, still ~1 min of sympy)
    basis, ref = _quadrature_check(L, alphas, 1 if L == 1 else 0)
    for op, M in ref.items():
        E = basis.mats[op]
        scale = np.abs(M).max()
        assert np.abs(E - M).max() < 1e-8 * scale, (L, op, np.abs(E - M).max(), scale)


@pytest.mark.parametrize("L", [1, 2, 3, 4])
def test_exchange_symmetry_union_equals_unsymmetrized(L):
    he = helium_like(2.0)
    bl = block(4, (1.3, 1.3, 2.1))
    eS = CoupledBasis(L, bl, "S").solve(he, 6, vectors=False)[0]
    eA = CoupledBasis(L, bl, "A").solve(he, 6, vectors=False)[0]
    eN = CoupledBasis(L, bl, None).solve(he, 6, vectors=False)[0]
    union = np.sort(np.concatenate([eS, eA]))[:6]
    assert np.allclose(union, eN, atol=1e-9), (union, eN)


def test_matrices_symmetric_and_mass_linear():
    fam = symmetric_family(7.0)
    b = CoupledBasis(2, block(4, (1.1, 1.1, 0.8)), "A")
    for m in b.mats.values():
        assert np.allclose(m, m.T, atol=1e-13 * np.abs(m).max())
    H = b.hamiltonian(fam)
    assert np.allclose(H, b.mats["T01"] / 7.0 + b.mats["T2"] + b.mats["V01"] - b.mats["V02+V12"])


@pytest.mark.parametrize("L", [3, 4])
def test_weighted_tensor_components_equal_full_traceless_contraction(L):
    """sum_K A_c,K A_d,K over the full 3^L traceless tensor (independent detracer formula, test code)
    equals the weighted multiset sum of angular.tensor_components (harmonic projection), for every c, d;
    and the full tensor built from tensor_components is traceless."""
    sub = {angular._a[k]: VA[k] for k in range(3)}
    sub.update({angular._b[k]: VB[k] for k in range(3)})
    comps = angular.components(L)
    for c in comps:
        Ac, wc = angular.tensor_components(c)
        for d in comps:
            Ad, _ = angular.tensor_components(d)
            lhs = sp.expand(sum(w * x * y for w, x, y in zip(wc, Ac, Ad)).subs(sub, simultaneous=True))
            Fc, Fd = _tensor_cart(c), _tensor_cart(d)
            rhs = sp.expand(sum(x * y for x, y in zip(Fc, Fd)))
            assert sp.expand(lhs - rhs) == 0, (c, d)


@pytest.mark.parametrize("L", [3])
def test_matrix_elements_against_quadrature_second_order_form_L3(L):
    alphas = (0.9, 1.4, 1.7)
    basis, ref = _quadrature_check(L, alphas, 0)
    for op, M in ref.items():
        E = basis.mats[op]
        scale = np.abs(M).max()
        assert np.abs(E - M).max() < 1e-8 * scale, (L, op, np.abs(E - M).max(), scale)
