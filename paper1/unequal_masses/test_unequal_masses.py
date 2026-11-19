"""Tests for paper1/unequal_masses (twochannel.py, massfunc_gu.py, run_exact.LeanBasis).  Written 2026-11-16.

    python3 -m pytest -p no:cacheprovider paper1/unequal_masses/test_unequal_masses.py -q   (~2-3 min)
"""
import os
import sys
import warnings

import numpy as np
import pytest
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")
os.environ.setdefault("THREEBODY_BLAS_THREADS", "4")

import twochannel as T  # noqa: E402


def root(model, lo, hi):
    return brentq(model.F, lo, hi, xtol=1e-11, rtol=1e-13)


# ---------------------------------------------------------------- r = 1 limits (symmetric family)
def test_r1_reproduces_coupling_py():
    """bo_offset/coupling.py: BO v=1 11.5814565935 (lam), BO+ad v=1 10.5631790317, BO+ad L=1 5.0705652790."""
    assert abs(root(T.Model("BO", 1.0, 0), 11.4, 11.7) - 11.5814565935) < 1e-8
    assert abs(root(T.Model("BOad", 1.0, 0), 10.4, 10.7) - 10.5631790317) < 1e-8
    assert abs(root(T.Model("BOad", 1.0, 1), 4.9, 5.2) - 5.0705652790) < 1e-8


def test_r1_mass_fold_matches_massfunc():
    """Non-perturbative symmetric offset 0.5281 (massfunc/analysis.txt, L=0 v=1)."""
    lam = root(T.Model("BOad+m", 1.0, 0), 9.9, 10.2)
    assert abs((10.5631790317 - lam) - 0.5281) < 2e-4


def test_bo_independent_of_r():
    assert abs(root(T.Model("BO", 2.0, 0), 11.4, 11.7) - 11.5814565935) < 1e-8


# ---------------------------------------------------------------- mass matrices
def test_mass_table_sum_rules_and_limits():
    tab = T.tables()
    df = tab.mass
    assert np.max(np.abs(df.mr_gg2 - 2.0)) < 1e-6           # (h - e_g)^-1 d_x g = -x g on complete Pi_u
    sel = df.R > 0.3
    assert np.max(np.abs(df.mv_gg2[sel] + df.mv_uu2[sel] - 4.0)) < 1e-6   # 2p sigma_u moved between rows
    m = tab.masses(39.9)
    assert abs(m["mv_gg0"] - 0.5) < 1e-3 and abs(m["mv_gu1"] + 2) < 2e-3 and abs(m["mr_gu1"] + 2) < 2e-3
    for r in (1.5, 2.0):
        b = T.UQ.beta_of(r)
        mm = np.array([[0.5 + 2 * b * b, -2 * b], [-2 * b, 0.5 + 2 * b * b]])
        assert abs(np.linalg.eigvalsh(mm)[0] - 2 / (1 + r) ** 2) < 1e-12   # open-channel reduced-mass value


def test_massfunc_gu_beta0_row_equals_symmetric():
    import massfunc_gu as MG
    q = MG.mass_matrices(6.0)
    assert abs(q["mv_gg0"] - 0.660321) < 2e-6 and abs(q["mr_gg0"] - 0.382848) < 2e-6


# ---------------------------------------------------------------- two-channel numerics
@pytest.mark.parametrize("kind", ["2ch", "2ch+m"])
def test_two_channel_numerical_parameters(kind):
    """lam_c(J=0, v=1, r=2) is insensitive to the matching radii and the inner start (to 1e-6)."""
    m = T.Model(kind, 2.0, 0)
    lo, hi = (16.2, 16.6) if kind == "2ch" else (15.9, 16.2)
    ref = root(m, lo, hi)
    old = (T.RX, m.rm)
    try:
        T.RX = 4.5
        m.rm = 1.6 * old[1]
        alt = root(m, lo, hi)
    finally:
        T.RX, m.rm = old
    assert abs(alt - ref) < 1e-6, (ref, alt)


def test_two_channel_small_beta_scaling():
    """Shift from BO+ad goes as sqrt(beta) for beta -> 0 (THEORY.md section 4)."""
    base = 10.5631790317
    s1 = root(T.Model("2ch", 1.001, 0), 10.6, 11.0) - base
    s2 = root(T.Model("2ch", 1.01, 0), 11.0, 11.5) - base
    b1, b2 = T.UQ.beta_of(1.001), T.UQ.beta_of(1.01)
    assert abs((s2 / s1) / np.sqrt(b2 / b1) - 1) < 0.1


# ---------------------------------------------------------------- exact-side lean basis
def test_lean_basis_matches_coupled_basis():
    sys.path.insert(0, os.path.join(HERE, "..", "..", "solvers", "threebody"))
    import run_exact as R
    from threebody.coupled import CoupledBasis
    sc = dict(M_opt=12.0, scales=[[4.0, 4.0, 1.3], [1.5, 1.5, 1.15]])
    t = "J0v1"
    bl = R.blocks_for(t, sc, 8, alpha_s=1e-2)
    lb = R.LeanBasis(0, bl, 2.0)
    cb = CoupledBasis(0, bl, None)
    for M in (8.0, 12.0):
        e1 = lb.energies(M, 3)
        e2 = cb.solve(R.system(M, 2.0), 3, vectors=False)[0]
        assert np.max(np.abs(e1 - e2)) < 1e-10
    assert R.threshold(10.0, 2.0) == -0.5 * 20 / 21
