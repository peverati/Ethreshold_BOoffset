"""Validation tests.  Run from the solvers/ folder:  python -m pytest bo_h2plus/tests -q

Reference values marked (memory) were quoted from memory by the requester and
are used with loose tolerances; the tight tests are internal cross-checks.
"""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from bo_h2plus import electronic as el
from bo_h2plus import potential as P
from bo_h2plus.adiabatic import D_of_R
from bo_h2plus.nuclear import solve_levels, count_levels, thresholds
from bo_h2plus.scan import critical_masses, HARTREE_CM, M_H2PLUS

# ---------------------------------------------------------------- electronic


def test_R2_reference():
    # (memory) -0.6026342144949 ; literature value -0.60263421449494
    U = el.total_energy(2.0)
    assert abs(U - (-0.6026342144949)) < 2e-12


@pytest.mark.parametrize("R", [0.05, 0.7, 2.0, 6.0, 25.0, 120.0])
def test_truncation_converged(R):
    e0 = el.electronic_energy(R)
    e1 = el.electronic_energy(R, n_ang=80 + int(2.5 * R), n_rad=140)
    assert abs(e0 - e1) < 2e-11


def _shoot_radial(R, nxi=40.0):
    """Independent radial check: integrate the xi ODE with solve_ivp for given p
    (lam from the Legendre angular problem) and find p where X stays bounded."""
    def tail(p):
        lam = el._lam_ang(p, 0, 60)
        a = 2 * R
        c = lam + a - p * p          # 2 X'(1) + c X(1) = 0 (regularity)
        x0 = 1e-7
        y0 = [1 - 0.5 * c * x0, -0.5 * c]

        def f(x, y):                 # x = xi - 1, (xi^2-1) = x (x+2)
            xi = 1 + x
            X, dX = y
            return [dX, -(2 * xi * dX + (lam + a * xi - p * p * xi * xi) * X) / (x * (x + 2))]
        s = solve_ivp(f, (x0, nxi / p), y0, method="DOP853", rtol=1e-13, atol=1e-300)
        # scale out the growing solution exp(+p xi)
        return s.y[0, -1] * math.exp(-p * (1 + nxi / p))
    pg = R * math.sqrt(-el.electronic_energy(R) / 2)
    p = brentq(tail, pg * (1 - 1e-4), pg * (1 + 1e-4), xtol=1e-14)
    return -2 * p * p / R ** 2


@pytest.mark.parametrize("R", [1.0, 2.0, 4.0])
def test_independent_radial_shooting(R):
    assert abs(_shoot_radial(R) - el.electronic_energy(R)) < 1e-9


def test_equilibrium():
    eq = el.equilibrium()
    assert abs(eq.Re - 1.9971933) < 1e-6          # (memory) 1.9972
    assert abs(eq.Emin - (-0.6026346191)) < 1e-9    # (memory) -0.6026346
    assert abs(el.dEdR(eq.Re)) < 1e-10


def test_united_atom_limit():
    for R in (1e-2, 1e-3, 1e-4):
        e = el.electronic_energy(R)
        # E_el = -2 + (8/3) R^2 + ...  (united-atom expansion)
        assert abs((e + 2) / R ** 2 - 8 / 3) < 0.05 + 30 * R


def test_separated_atom_polarization():
    # (U + 1/2) R^4 -> -9/4 with next term -15/(2 R^2) ; fit both
    Rs = np.array([60.0, 80.0, 100.0, 140.0, 200.0])
    y = np.array([(el.total_energy(R) + 0.5) * R ** 4 for R in Rs])
    A = np.vstack([np.ones_like(Rs), Rs ** -2, Rs ** -3]).T
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    assert abs(c[0] + 2.25) < 1e-4
    assert abs(c[1] + 7.5) < 0.3


@pytest.mark.parametrize("R", [0.8, 2.0, 5.0])
def test_virial_euler(R):
    """E_el(R,Z) = Z^2 f(ZR)  =>  2 E_el - Z dE/dZ + R dE/dR = 0 (virial theorem)."""
    h = 1e-4
    E = el.electronic_energy(R)
    dZ = (el.electronic_energy(R, Z=1 + h) - el.electronic_energy(R, Z=1 - h)) / (2 * h)
    dR = (el.electronic_energy(R + h * R) - el.electronic_energy(R - h * R)) / (2 * h * R)
    assert abs(2 * E - dZ + R * dR) < 1e-7


@pytest.mark.parametrize("mu", [2 / 3, 0.9, 0.9997])
def test_mass_scaling(mu):
    """Electron mass mu  <=>  charge Z = mu, energy/mu:  E(R;mu) = E(R;Z=mu)/mu = mu E(mu R;1)."""
    R = 2.0
    a = el.electronic_energy(R, Z=mu) / mu
    b = el.electronic_energy_mu(R, mu)
    assert abs(a - b) < 1e-12


def test_2psigma_u():
    # (memory) 2p sigma_u total energy at R = 2: -0.16753439
    assert abs(el.total_energy(2.0, "2pu") - (-0.16753439)) < 1e-7


def test_table_interpolation():
    rng = np.random.default_rng(0)
    for R in np.exp(rng.uniform(np.log(0.011), np.log(199.0), 25)):
        assert abs(P.E_el_interp(R)[0] - el.electronic_energy(R)) < 2e-11


# ---------------------------------------------------------------- adiabatic


def test_adiabatic_internal():
    d = D_of_R(2.0, details=True)
    assert abs(d["norm"] - 1) < 1e-12 and abs(d["overlap"]) < 1e-10
    assert abs(D_of_R(2.0, na=80, nb=80) - d["D"]) < 1e-12
    assert abs(D_of_R(2.0, method="cartesian", h=1e-3) - d["D"]) < 1e-8


def test_adiabatic_large_R():
    d = D_of_R(40.0, details=True)
    assert abs(d["D"] - 0.25) < 2e-6
    assert abs(d["radial"] - 1 / 12) < 1e-4 and abs(d["angular"] - 1 / 6) < 1e-4


# ---------------------------------------------------------------- nuclear


@pytest.fixture(scope="module")
def h2p():
    return solve_levels(M_H2PLUS, 0)


def test_h2plus_levels(h2p):
    lv = h2p.levels
    assert len(lv) == 20 == count_levels(M_H2PLUS, 0)           # (memory) v = 0..19
    assert max(l.err_E for l in lv) < 1e-9
    fund = (lv[1].E - lv[0].E) * HARTREE_CM
    assert abs(fund - 2191.1) < 0.5                            # (memory) ~2191 cm-1
    assert -0.5975 < lv[0].E < -0.5970


def test_h2plus_vs_nonBO(h2p):
    exact = -0.597139063                                       # (memory) non-BO ground state
    E_bo = h2p.levels[0].E
    E_ad = solve_levels(M_H2PLUS, 0, adiabatic=True, check=False).levels[0].E
    assert -2e-4 < E_bo - exact < 0       # BO lies below
    assert abs(E_ad - exact) < 2e-6       # BO + adiabatic within the non-adiabatic correction


def test_h2plus_rotation_and_harmonic(h2p):
    s1 = solve_levels(M_H2PLUS, 1, check=False)
    d = (s1.levels[0].E - h2p.levels[0].E) * HARTREE_CM
    assert abs(d / 2 - 29.1) < 0.3                             # B_0 (memory: ~30)
    eq = el.equilibrium(mu=P.mu_e(M_H2PLUS))
    om = math.sqrt(eq.k / (M_H2PLUS / 2)) * HARTREE_CM
    assert abs(om - 2322) < 3                                  # (memory) ~2322


def test_grid_convergence(h2p):
    fine = solve_levels(M_H2PLUS, 0, ppw=3.0, box=2.0, check=False)
    for a, b in zip(h2p.levels, fine.levels):
        assert abs(a.E - b.E) < 1e-10
        assert abs(a.R_mean - b.R_mean) < 1e-6


def test_density_normalized(h2p):
    from numpy import trapz
    for v in (0, 5, 19):
        p = h2p.P[v]
        assert abs(trapz(p, h2p.R) - 1) < 1e-3   # trapezoid on a nonuniform grid


@pytest.mark.parametrize("M", [1.0, 5.0, 50.0, 1e4])
def test_count_vs_dvr(M):
    assert solve_levels(M, 0, check=False).n_dvr == count_levels(M, 0)


def test_V_BO_mass_scaling():
    """Interpolated V_BO(R; M) = mu_e U(mu_e R) agrees with the direct solver at mass mu_e."""
    M = 3.0
    me = P.mu_e(M)
    for R in (1.0, 2.5, 7.0):
        direct = el.electronic_energy_mu(R, me) + 1 / R
        assert abs(P.V_BO(R, M)[0] - direct) < 2e-11


def test_thresholds():
    bo, ex = thresholds(M_H2PLUS)
    assert abs(ex - bo - M_H2PLUS / ((2 * M_H2PLUS + 1) * (2 * M_H2PLUS + 2))) < 1e-15


def test_critical_masses():
    res = critical_masses(3)
    for d in res:
        assert d["est_abs_err"] < 1e-7 * d["M_c"]
        n, Mc = d["n"], d["M_c"]
        below = solve_levels(Mc * (1 - 2e-3), 0, check=False).n_dvr
        above = solve_levels(Mc * (1 + 2e-3), 0, check=False).n_dvr
        assert below == n - 1 and above == n
