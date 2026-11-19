"""Verification of the kinetic-energy operator and of the exact matrix elements.

1. The distance-coordinate Laplacian used in the derivation is checked against a
   brute-force Cartesian Laplacian (sympy, 9 Cartesian coordinates) for a generic
   test function f(r01, r02, r12) at random configurations.
2. Matrix elements produced by the engine (symmetric, first-derivative form) are
   checked against direct numerical quadrature of <F_a | -(1/2m_i) Lap_i F_b>
   (second-order form) and of the overlap/Coulomb integrals, for unequal masses
   and unequal scale parameters.
"""
import numpy as np
import sympy as sp
from scipy.special import roots_laguerre

from threebody import PerimetricBasis
from threebody.operators import kinetic_terms, overlap_terms, power_terms, PAIRS

R01, R02, R12 = sp.symbols("R01 R02 R12", positive=True)
DIST = {(0, 1): R01, (0, 2): R02, (1, 2): R12}


def _test_function():
    # generic, no symmetry, contains all three distances and cross terms
    return (1 + R01 * R02 / 3 + R12**2 / 5 + R01 * R12 / 7) * sp.exp(-0.7 * R01 - 1.1 * R02 - 0.3 * R12)


def _laplacian_distance_form(f, i):
    j, k = [x for x in (0, 1, 2) if x != i]
    rij, rik, rjk = DIST[tuple(sorted((i, j)))], DIST[tuple(sorted((i, k)))], DIST[tuple(sorted((j, k)))]
    cos_i = (rij**2 + rik**2 - rjk**2) / (2 * rij * rik)
    return (sp.diff(f, rij, 2) + 2 / rij * sp.diff(f, rij) + sp.diff(f, rik, 2) + 2 / rik * sp.diff(f, rik)
            + 2 * cos_i * sp.diff(f, rij, rik))


def test_laplacian_against_cartesian():
    X = sp.symbols("x0:9", real=True)
    pos = [sp.Matrix(X[3 * n:3 * n + 3]) for n in range(3)]
    d = {p: sp.sqrt(((pos[p[0]] - pos[p[1]]).T * (pos[p[0]] - pos[p[1]]))[0]) for p in PAIRS}
    f = _test_function()
    fc = f.subs({R01: d[(0, 1)], R02: d[(0, 2)], R12: d[(1, 2)]})
    rng = np.random.default_rng(1)
    for i in range(3):
        lap_cart = sum(sp.diff(fc, X[3 * i + c], 2) for c in range(3))
        lap_dist = _laplacian_distance_form(f, i).subs({R01: d[(0, 1)], R02: d[(0, 2)], R12: d[(1, 2)]})
        fa = sp.lambdify(X, lap_cart, "numpy")
        fb = sp.lambdify(X, lap_dist, "numpy")
        for _ in range(5):
            x = rng.normal(size=9)
            a, b = fa(*x), fb(*x)
            assert abs(a - b) < 1e-11 * max(1.0, abs(a)), (i, a, b)


def _quadrature_matrix(basis, operator_on_ket, measure):
    """Brute-force <F_a|O F_b> on a tensor Gauss-Laguerre grid in perimetric coordinates."""
    s = sp.symbols("s0:3", nonnegative=True)
    subs = {R01: (s[0] + s[1]) / 2, R02: (s[0] + s[2]) / 2, R12: (s[1] + s[2]) / 2}
    nq = 24
    x, w = roots_laguerre(nq)
    grids = []
    for m in range(3):
        al = basis.alphas[m]
        grids.append((x / al, w * np.exp(x) / al))  # nodes, weights for int ds (exp is in the f's)
    S0, S1, S2 = np.meshgrid(grids[0][0], grids[1][0], grids[2][0], indexing="ij")
    W = np.einsum("i,j,k->ijk", grids[0][1], grids[1][1], grids[2][1])
    funcs = []
    for a in basis.index:
        expr = 1
        for m in range(3):
            expr *= sp.assoc_laguerre(int(a[m]), 0, basis.alphas[m] * s[m]) * sp.exp(-basis.alphas[m] * s[m] / 2)
        funcs.append(expr)
    vals = [sp.lambdify(s, f, "numpy")(S0, S1, S2) for f in funcs]
    meas = sp.lambdify(s, measure.subs(subs), "numpy")(S0, S1, S2)
    ovals = []
    for f in funcs:
        fr = f.subs({s[0]: R01 + R02 - R12, s[1]: R01 + R12 - R02, s[2]: R02 + R12 - R01})
        o = operator_on_ket(fr).subs(subs)
        ovals.append(sp.lambdify(s, o, "numpy")(S0, S1, S2) * np.ones_like(S0))
    n = len(funcs)
    M = np.empty((n, n))
    for p in range(n):
        for q in range(n):
            M[p, q] = np.sum(W * vals[p] * ovals[q] * meas)
    return M


def test_matrix_elements_against_quadrature():
    basis = PerimetricBasis(2, (0.9, 1.4, 1.7), None)
    w = R01 * R02 * R12 / 4  # measure incl. Jacobian 1/4
    S_ref = _quadrature_matrix(basis, lambda f: f, w)
    assert np.allclose(basis.matrix(overlap_terms()), S_ref, rtol=1e-10, atol=1e-12)
    for (i, j) in PAIRS:
        V_ref = _quadrature_matrix(basis, lambda f, i=i, j=j: f / DIST[(i, j)], w)
        assert np.allclose(basis.matrix(power_terms(i, j, -1)), V_ref, rtol=1e-10, atol=1e-12)
    for i in range(3):
        T_ref = _quadrature_matrix(basis, lambda f, i=i: -sp.Rational(1, 2) * _laplacian_distance_form(f, i), w)
        assert np.allclose(basis.matrix(kinetic_terms(i)), T_ref, rtol=1e-9, atol=1e-11), i
