"""Tests of the multi-scale (multi-block) basis in threebody/threshold.py."""
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh

from threebody import PerimetricBasis, solve, symmetric_family, ThreeBodySystem
from threebody import laguerre
from threebody.threshold import (Block, MultiBlockBasis, compact_block, diffuse_block, one_d_mixed,
                                 raw_matrix, family_operators, total_degree_index)


def _f(n, a, s, deriv=False):
    """f_n(s;a) and derivative by brute force (numpy Laguerre)."""
    from numpy.polynomial import laguerre as L
    c = np.zeros(n + 1); c[n] = 1
    x = a * s
    v = L.lagval(x, c) * np.exp(-x / 2)
    if not deriv:
        return v
    return a * (L.lagval(x, L.lagder(c)) - 0.5 * L.lagval(x, c)) * np.exp(-x / 2)


def test_mixed_one_d_against_quad():
    a, b = 1.3, 0.07
    for p in (0, 1, 3):
        for da in (0, 1):
            for db in (0, 1):
                M = one_d_mixed(a, 5, b, 7, p, da, db)
                for (i, j) in [(0, 0), (2, 5), (5, 7), (3, 1)]:
                    ref = quad(lambda s: _f(i, a, s, da) * _f(j, b, s, db) * s ** p, 0, np.inf, limit=400)[0]
                    assert abs(M[i, j] - ref) < 1e-9 * max(1, abs(ref))


def test_mixed_continuous_to_same_scale():
    a = 0.8
    for p in (0, 2):
        for da, db in [(0, 0), (1, 0), (1, 1)]:
            ex = laguerre.one_d_matrix(6, p, da, db, a)
            mx = one_d_mixed(a, 6, a * (1 + 1e-9), 6, p, da, db)
            assert np.max(np.abs(ex - mx)) < 1e-6 * np.max(np.abs(ex))


def test_compact_block_reproduces_package():
    M, N, al = 10.0, 10, (0.9, 0.9, 1.1)
    s = symmetric_family(M)
    ref = solve(s, PerimetricBasis(N, al, "S"), nstates=3).energies
    mb = MultiBlockBasis([compact_block(N, al[0], al[2])], ops=("S", "TX", "Te", "V01", "VXe"))
    e = mb.energies(M, 3, tol=1e-13)
    assert np.allclose(e, ref, atol=1e-11, rtol=0)


def test_symmetrization_against_unsymmetrized_union():
    """phi = F +/- PF construction vs. raw diagonalization of block U mirror (no symmetry used)."""
    M = 7.0
    B = Block((1.0, 0.3, 0.9), total_degree_index(5))
    ops = family_operators()
    blocks = [B, B.mirror()]
    def full(name):
        cache = {}
        rows = [np.hstack([raw_matrix(ops[name], X, Y, cache) for Y in blocks]) for X in blocks]
        return np.vstack(rows)
    S = full("S"); H = full("TX") / M + full("Te") + full("V01") - full("VXe")
    e_all = eigh(H, S, eigvals_only=True)
    eS = MultiBlockBasis([B], ops=("S", "TX", "Te", "V01", "VXe"), sign=+1).energies(M, 4, tol=1e-14)
    eA = MultiBlockBasis([B], ops=("S", "TX", "Te", "V01", "VXe"), sign=-1).energies(M, 4, tol=1e-14)
    merged = np.sort(np.concatenate([eS, eA]))[:4]
    assert np.allclose(merged, e_all[:4], atol=1e-9)
    assert eS[0] < eA[0]


def test_union_is_variational_and_hf_derivative():
    M = 12.0
    C = compact_block(8, 0.8, 1.0)
    D = diffuse_block(0.92, 0.1, 2, 8)
    e1 = MultiBlockBasis([C]).energies(M, 2)
    mb = MultiBlockBasis([C, D])
    e2 = mb.energies(M, 2)
    assert np.all(e2 <= e1 + 1e-12)
    ob = mb.observables(M, state=1)
    h = 1e-4
    fd = (mb.energies(M + h, 2)[1] - mb.energies(M - h, 2)[1]) / (2 * h)
    assert abs(fd - ob["dE_dM"]) < 1e-7
