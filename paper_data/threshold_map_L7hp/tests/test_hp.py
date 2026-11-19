"""Tests of the double-double assembly (hp_assembly.py).  Written 2026-11-17.
Run from paper_data/threshold_map_L7hp:  python -m pytest -q tests"""
import os
import sys
from fractions import Fraction

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import hp_assembly as H  # noqa: E402
from threebody import laguerre  # noqa: E402
from threebody.coupled import CoupledBasis, IdxBlock, block  # noqa: E402
from threebody.system import ThreeBodySystem, symmetric_family  # noqa: E402
from threebody.threshold import diffuse_index  # noqa: E402

rng = np.random.default_rng(1)


def exact(h, l):
    return Fraction(float(h)) + Fraction(float(l))


def test_dd_primitives_exact():
    a, b = rng.normal(size=50) * 1e3, rng.normal(size=50) * 1e-3
    al, bl = a * 1e-17, b * 1e-17
    sh, sl = H.dd_add(a, al, b, bl)
    ph, pl = H.dd_mul(a, al, b, bl)
    for i in range(50):
        xs = exact(a[i], al[i]) + exact(b[i], bl[i])
        xp = exact(a[i], al[i]) * exact(b[i], bl[i])
        assert abs(exact(sh[i], sl[i]) - xs) <= abs(xs) * Fraction(1, 10 ** 30)
        assert abs(exact(ph[i], pl[i]) - xp) <= abs(xp) * Fraction(1, 10 ** 29)


@pytest.mark.parametrize("p,da,db", [(0, 0, 0), (3, 1, 0), (7, 0, 1), (12, 1, 1), (16, 0, 0)])
def test_same_scale_matches_double(p, da, db):
    al, n = 1.7, 25
    hi, lo = H.one_d_dd(al, n, al, n, p, da, db)
    ref = laguerre.one_d_matrix(n, p, da, db, al)
    assert np.max(np.abs(hi - ref)) <= 1e-13 * np.max(np.abs(ref))   # double itself is inexact once |J^p| > 2^53


def test_same_scale_exact_reference():
    """Q (exact integers) against an independent Fraction computation of alpha^-(p+1) J^p with derivatives."""
    n, p, al = 10, 9, Fraction(3, 2)
    size = n + 1 + p
    J = [[Fraction(0)] * size for _ in range(size)]
    for i in range(size):
        J[i][i] = Fraction(2 * i + 1)
        if i + 1 < size:
            J[i][i + 1] = J[i + 1][i] = Fraction(-(i + 1))
    M = [[Fraction(int(i == j)) for j in range(size)] for i in range(size)]
    for _ in range(p):
        M = [[sum(M[i][k] * J[k][j] for k in range(size)) for j in range(size)] for i in range(size)]
    S = [[M[i][j] / al ** (p + 1) for j in range(n + 1)] for i in range(n + 1)]
    D = [[Fraction(-1) if k < m else (Fraction(-1, 2) if k == m else Fraction(0)) for m in range(n + 1)] for k in range(n + 1)]
    Sd = [[al * sum(D[k][i] * S[k][j] for k in range(n + 1)) for j in range(n + 1)] for i in range(n + 1)]
    hi, lo = H.one_d_dd(float(al), n, float(al), n, p, 1, 0)
    for i in range(n + 1):
        for j in range(n + 1):
            assert abs(exact(hi[i, j], lo[i, j]) - Sd[i][j]) <= abs(Sd[i][j]) * Fraction(1, 10 ** 30) + Fraction(1, 10 ** 40)


def test_mixed_scale_quadrature_vs_expansion():
    for (na, a, nb, b, p, da, db) in [(6, 0.7, 9, 2.3, 0, 0, 0), (8, 1.1, 5, 0.03, 3, 1, 0), (7, 0.4, 7, 5.0, 5, 1, 1)]:
        qh, ql = H.mixed_moment_dd(na, a, nb, b, p, da, db)
        eh, el = H.mixed_moment_exact_mp(na, a, nb, b, p, da, db)
        q = np.vectorize(exact)(qh, ql)
        e = np.vectorize(exact)(eh, el)
        scale = max(abs(x) for x in e.ravel())
        assert max(abs(x) for x in (q - e).ravel()) <= scale * Fraction(1, 10 ** 28)


def test_mixed_scale_matches_double():
    hi, _ = H.one_d_dd(0.6, 12, 1.9, 30, 4, 1, 0)
    S = laguerre.mixed_moment_matrix(12, 0.6, 30, 1.9, 4)
    ref = 0.6 * (laguerre.derivative_matrix(13).T @ S)
    assert np.max(np.abs(hi - ref)) <= 1e-12 * np.max(np.abs(ref))


@pytest.mark.parametrize("L,sym", [(1, "A"), (2, "S")])
def test_coupled_matches_double_low_L(L, sym):
    bl = block(6, (1.3, 1.3, 1.1))
    bd = CoupledBasis(L, bl, sym)
    bh = H.HPCoupledBasis(L, bl, sym)
    for op in bd.mats:
        assert np.max(np.abs(bd.mats[op] - bh.mats[op])) <= 1e-13 * np.max(np.abs(bh.mats[op]))
    s = symmetric_family(5.0)
    ed = bd.solve(s, 3, 1e-10, vectors=False)[0]
    eh = bh.solve(s, 3, 1e-10, vectors=False)[0]
    assert np.max(np.abs(ed - eh)) < 1e-10


def test_mixed_blocks_unsymmetrized_matches_double():
    blocks = [block(4, (1.4, 1.2, 1.0)), IdxBlock((0.9, 0.05, 0.9), diffuse_index(2, 10)[:, [1, 0, 2]], "O")]
    bd = CoupledBasis(1, blocks, None)
    bh = H.HPCoupledBasis(1, blocks, None)
    for op in bd.mats:
        assert np.max(np.abs(bd.mats[op] - bh.mats[op])) <= 1e-12 * np.max(np.abs(bh.mats[op]))
    s = ThreeBodySystem((7.0, 14.0, 1.0), (1.0, 1.0, -1.0))
    ed = bd.solve(s, 2, 1e-10, vectors=False)[0]
    eh = bh.solve(s, 2, 1e-10, vectors=False)[0]
    assert np.max(np.abs(ed - eh)) < 1e-9


def test_parallel_equals_serial(monkeypatch):
    bl = block(5, (1.3, 1.3, 1.1))
    monkeypatch.setattr(H, "NPROC", 1)
    b1 = H.HPCoupledBasis(2, bl, "S")
    monkeypatch.setattr(H, "NPROC", 4)
    b4 = H.HPCoupledBasis(2, bl, "S")
    for op in b1.mats:
        assert np.array_equal(b1.mats[op], b4.mats[op])


def test_two_pass_projection_orthonormal(monkeypatch):
    monkeypatch.setattr(H.HPCoupledBasis, "two_pass", True)
    b = H.HPCoupledBasis(2, block(8, (1.3, 1.3, 1.1)), "S")
    pr = b.project(1e-14)
    assert pr["ortho_err"] < 1e-8           # ~eps / tol amplification; guard threshold
    e2 = b.solve(symmetric_family(5.0), 2, 1e-14, vectors=False)[0]
    b1 = H.HPCoupledBasis(2, block(8, (1.3, 1.3, 1.1)), "S")
    b1.two_pass = False
    e1 = b1.solve(symmetric_family(5.0), 2, 1e-10, vectors=False)[0]
    assert np.all(e2 <= e1 + 1e-12)        # larger kept space, still variational
