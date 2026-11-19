"""Energies against literature targets and exact limits (moderate basis sizes, < ~1 min total)."""
import numpy as np
import pytest

from threebody import (ThreeBodySystem, PerimetricBasis, solve, helium_like, symmetric_family,
                       family_threshold, M_PROTON)

# scales optimized at N = 10 (see README); fixed here so the tests are fast
HE_ALPHAS = (2.6562541, 2.6562541, 4.52188892)
HM_ALPHAS = (0.74411999, 0.74411999, 1.61470616)
PS_ALPHAS = (0.36950233, 0.36950233, 0.75592905)
H2P_ALPHAS = (6.69517288, 6.69517288, 1.30638368)


@pytest.mark.parametrize("system,alphas,ref,tol", [
    (helium_like(2.0), HE_ALPHAS, -2.903724377034, 1e-9),
    (helium_like(1.0), HM_ALPHAS, -0.527751016544, 1e-8),
    (symmetric_family(1.0), PS_ALPHAS, -0.262005070233, 1e-8),
])
def test_literature_ground_states(system, alphas, ref, tol):
    e = solve(system, PerimetricBasis(20, alphas, "S"), 1).energies[0]
    assert e >= ref - 1e-11  # variational: never below the (accurate) reference
    assert e - ref < tol


def test_h2plus_ground_state_moderate():
    s = symmetric_family(M_PROTON)
    e = solve(s, PerimetricBasis(20, H2P_ALPHAS, "S"), 2).energies
    assert -0.597139063 - 1e-9 < e[0] < -0.597139063 + 2e-6
    # second root (v=1) is far below threshold: rigorous evidence of >= 2 bound S states
    assert e[1] < family_threshold(M_PROTON) - 0.05


def test_variational_monotonic():
    s = helium_like(1.0)
    es = [solve(s, PerimetricBasis(N, HM_ALPHAS, "S"), 3).energies for N in (6, 8, 10, 12)]
    for a, b in zip(es, es[1:]):
        assert np.all(b <= a + 1e-13)


def test_separable_limit_exact():
    # infinite nucleus, electron-electron interaction switched off: hydrogenic products
    Z = 2.0
    s = ThreeBodySystem((1, 1, np.inf), (-1, -1, Z), couplings={(0, 1): 0.0}, identical=(0, 1))
    e = solve(s, PerimetricBasis(3, (Z, Z, 2 * Z), "S"), 1).energies[0]
    assert abs(e - (-Z * Z)) < 1e-12
    # unequal "electron" masses, no symmetry: E = -Z^2 (m0 + m1) / 2
    m0, m1 = 1.0, 3.0
    s = ThreeBodySystem((m0, m1, np.inf), (-1, -1, Z), couplings={(0, 1): 0.0})
    e = solve(s, PerimetricBasis(3, (Z * m0, Z * m1, Z * (m0 + m1)), None), 1).energies[0]
    assert abs(e - (-Z * Z * (m0 + m1) / 2)) < 1e-12
    # 1s2s singlet-like excited state is approached from above
    s = ThreeBodySystem((1, 1, np.inf), (-1, -1, Z), couplings={(0, 1): 0.0}, identical=(0, 1))
    e = solve(s, PerimetricBasis(14, (Z, Z, 2 * Z), "S"), 2).energies
    assert -2.5 < e[1] < -2.5 + 1e-3


def test_mass_scaling_identity():
    # scaling all masses by lam scales energies by lam (with scales alpha scaled by lam)
    lam, M = 2.5, 7.0
    b1 = PerimetricBasis(10, (1.0, 1.0, 1.3), "S")
    b2 = PerimetricBasis(10, (lam, lam, 1.3 * lam), "S")
    e1 = solve(ThreeBodySystem((M, M, 1), (1, 1, -1), identical=(0, 1)), b1, 3).energies
    e2 = solve(ThreeBodySystem((lam * M, lam * M, lam), (1, 1, -1), identical=(0, 1)), b2, 3).energies
    assert np.allclose(e2, lam * e1, rtol=1e-11, atol=1e-13)


def test_relabelling_invariance_unequal_masses():
    # {M, M, 1} solved with the symmetry-adapted code (particles 0,1 = X) and with the
    # general unsymmetrized code with the electron relabelled as particle 1.
    M = 5.0
    al = (1.3, 1.3, 1.1)
    N = 8
    es = solve(symmetric_family(M), PerimetricBasis(N, al, "S"), 1).energies[0]
    ea = solve(symmetric_family(M), PerimetricBasis(N, al, "A"), 1).energies[0]
    # relabel: new 0 = X, new 1 = e, new 2 = X ; perimetric coords follow the vertices
    s = ThreeBodySystem((M, 1.0, M), (1, -1, 1))
    e_full = solve(s, PerimetricBasis(N, (al[0], al[2], al[1]), None), 2).energies
    assert abs(e_full[0] - min(es, ea)) < 1e-11
    assert es < ea


def test_threshold_and_bound_count_ps_minus():
    s = symmetric_family(1.0)
    r = solve(s, PerimetricBasis(16, PS_ALPHAS, "S"), 3)
    assert r.threshold == -0.25
    assert r.n_below_threshold() == 1
