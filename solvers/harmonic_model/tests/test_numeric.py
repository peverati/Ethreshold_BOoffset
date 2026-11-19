"""Closed forms versus independent numerics (FD diagonalization, quadrature, SVD)."""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

from harmonic_model import model_a as A, model_b as B, twobody as T, numerics as N, common as C

MASSES = [0.1, 1.0, 10.0, 1836.0]


# ------------------------------------------------------------------ energy
@pytest.mark.parametrize("M", MASSES)
@pytest.mark.parametrize("g", [0.0, 1.0])
def test_energy_and_hellmann_feynman(M, g):
    en = N.expectations_numeric(M, g)
    assert abs(en["E"] - B.energy(M, g)) / B.energy(M, g) < 2e-5
    # Hellmann-Feynman: closed-form derivative = - <T_R> - <T_rho>/(2M+1) (numerical T)
    hf_num = -en["T_R"] - en["T_rho"] / (2 * M + 1)
    assert abs(hf_num - B.dE_dlnM(M, g)) < 2e-4 * abs(B.dE_dlnM(M, g)) + 1e-6
    assert abs(B.hellmann_feynman(M, g) - B.dE_dlnM(M, g)) < 1e-12


@pytest.mark.parametrize("M", MASSES)
def test_model_a_closed_forms(M):
    k, kp = 1.3, -0.4
    assert abs(A.energy(M, k, kp) - B.energy(M, 0.0, k, kp)) < 1e-12
    assert abs(A.dE_dlnM(M, k, kp) - A.hellmann_feynman(M, k, kp)) < 1e-12
    h = 1e-3
    fd2 = (A.energy(M * np.exp(h), k, kp) - 2 * A.energy(M, k, kp) + A.energy(M * np.exp(-h), k, kp)) / h ** 2
    assert abs(fd2 - A.d2E_dlnM2(M, k, kp)) < 1e-5 * (1 + abs(fd2))
    w = N.normal_mode_frequencies(M, k, kp)
    p = A.params(M, k, kp)
    assert np.allclose(np.sort(w), np.sort([0, p["omR"], p["omrho"]]), atol=1e-6)
    with pytest.raises(ValueError):
        A.params(M, 1.0, -0.5)        # k' = -k/2: unbound


def test_rotational_tower():
    M, g = 50.0, 1.0
    for J in (1, 2):
        gr = N.RadialGrid(5.0, 8000)
        E = gr.ground(1 / M, lambda r: r ** 2 / 4 + (g + J * (J + 1) / M) / r ** 2)[0]
        exact = B.levels(M, 0, J, g) - 1.5 * B.params(M, g)["omrho"]
        assert abs(E - exact) < 1e-5


# --------------------------------------------------------------- moments
@pytest.mark.parametrize("M", MASSES)
def test_pair_moments(M):
    en = N.expectations_numeric(M, 1.0)
    pm = B.pair_moments(M, 1.0)
    assert abs(en["R"] - pm["R"]) / pm["R"] < 3e-5  # FD error near R=0 when lam is small
    assert abs(en["R2"] - pm["R2"]) / pm["R2"] < 5e-5
    pa = A.pair_moments(M)
    assert abs(pa["eta"] - A.LINDEMANN) < 1e-12
    assert abs(B.pair_moments(M, 0.0)["eta"] - A.LINDEMANN) < 1e-12
    I = quad(lambda r: B.pair_distribution(r, M, 1.0), 0, 40 / np.sqrt(B.params(M, 1.0)["a"]))[0]
    assert abs(I - 1) < 1e-8


# -------------------------------------------------------------- densities
@pytest.mark.parametrize("M", [0.3, 1.0, 10.0])
def test_cm_density(M):
    p = B.params(M, 1.0)
    y = np.linspace(1e-6, 12 / np.sqrt(p["a"]), 6000)
    nm = B.heavy_density_mid(y, M)
    s2 = 1 / (2 * p["b"] * (2 * M + 1) ** 2)
    xs = np.array([0.0, 0.3, 0.7])
    assert np.allclose(B.heavy_density_cm(xs, M), N.cm_density_numeric(xs, y, nm, s2), rtol=1e-6)
    # g = 0 reduces to the Model A Gaussian
    assert np.allclose(B.heavy_density_cm(xs, M, 0.0), A.heavy_density(xs, M), rtol=1e-10)


def test_critical_mass_cm():
    Mc = B.critical_mass_cm(1.0)
    for M, sign in ((Mc * 0.99, -1), (Mc * 1.01, 1)):
        p = B.params(M, 1.0)
        y = np.linspace(1e-7, 14 / np.sqrt(p["a"]), 20000)
        nm = B.heavy_density_mid(y, M)
        s2 = 1 / (2 * p["b"] * (2 * M + 1) ** 2)
        hh = 0.02 * np.sqrt(s2)
        n = N.cm_density_numeric([0, hh], y, nm, s2)
        curv = (n[1] - n[0]) / hh ** 2 / n[0]
        assert np.sign(curv) == sign
        assert abs(curv - B.cm_curvature_coefficient(M, 1.0)) < 2e-3 * abs(1 / s2)
    assert abs(B.cm_curvature_coefficient(Mc, 1.0)) < 1e-9


# ---------------------------------------------------------------- fidelity
def test_overlap_formulas_by_quadrature():
    for lam1, a1, lam2, a2 in ((0.3, 0.7, 2.1, 1.9), (0.0, 1.0, 0.0, 3.0)):
        num = quad(lambda r: 4 * np.pi * r ** 2 * C.chi_value(lam1, a1, r) * C.chi_value(lam2, a2, r), 0, 40)[0]
        assert abs(num - C.chi_overlap(lam1, a1, lam2, a2)) < 1e-10
    b1, b2 = 0.6, 2.2
    num = quad(lambda r: 4 * np.pi * r ** 2 * (b1 ** 3 / np.pi) ** .5 * np.exp(-b1 * r)
               * (b2 ** 3 / np.pi) ** .5 * np.exp(-b2 * r), 0, 80)[0]
    assert abs(num - T.coulomb_overlap(b1, b2)) < 1e-10


@pytest.mark.parametrize("M", [0.1, 1.0, 10.0])
def test_fidelity_and_inv_vs_fd_states(M):
    M2 = 1.7 * M
    for g in (0.0, 1.0):
        assert abs(N.overlap_numeric(M, M2, g) - B.fidelity(M, M2, g)) < 2e-5
        # dilation-invariant: maximize numerical overlap over s independently
        r = minimize_scalar(lambda t: -N.overlap_numeric(M, M2, g, s=np.exp(t), n=3000),
                            bounds=(-0.6, 0.6), method="bounded", options=dict(xatol=1e-6))
        assert abs(-r.fun - B.fidelity_inv(M, M2, g)) < 5e-5
    assert abs(B.fidelity_inv(M, M2, 0.0) - A.fidelity_inv(M, M2)) < 1e-10
    assert abs(B.fidelity(M, M2, 0.0) - A.fidelity(M, M2)) < 1e-12
    assert abs(B.fidelity_sep(M, M2, 0.0) - 1) < 1e-10


@pytest.mark.parametrize("M", MASSES)
@pytest.mark.parametrize("g", [0.0, 1.0])
def test_susceptibilities_vs_fd(M, g):
    cn = N.chi_numeric(M, g)
    cc = B.chi_parts(M, g)
    for key in ("chi_fixed", "chi_inv", "chi_sep"):
        assert abs(cn[key] - cc[key]) < 2e-5 * (1 + 10 * cc[key]), key
    # susceptibility = curvature of the closed-form fidelity
    d = 1e-3
    F = B.fidelity(M * np.exp(-d), M * np.exp(d), g)
    assert abs(2 * (1 - F) / (2 * d) ** 2 - cc["chi_fixed"]) < 1e-5
    Fi = B.fidelity_inv(M * np.exp(-0.01), M * np.exp(0.01), g)
    assert abs(2 * (1 - Fi) / 0.02 ** 2 - cc["chi_inv"]) < 2e-4
    if g == 0.0:
        assert abs(cc["chi_fixed"] - A.chi_fixed(M)) < 1e-12
        assert abs(cc["chi_inv"] - A.chi_inv(M)) < 1e-12
        assert abs(A.chi_fixed(M) - A.chi_inv(M) - A.chi_removed(M)) < 1e-14


def test_limits():
    assert abs(A.chi_fixed(1e-12) - 3 / 16) < 1e-10 and abs(A.chi_inv(1e-12)) < 1e-10
    assert abs(A.chi_fixed(1.0) - 5 / 48) < 1e-14 and abs(A.chi_inv(1.0) - 1 / 48) < 1e-14
    assert abs(A.chi_fixed(1e12) - 3 / 32) < 1e-10 and abs(A.chi_inv(1e12) - 3 / 64) < 1e-10
    c = B.chi_parts(1e12, 1.0)
    for key in ("chi_fixed", "chi_inv", "chi_sep"):
        assert abs(c[key] - 1 / 32) < 1e-5
    c = B.chi_parts(1e-10, 1.0)
    assert abs(c["chi_fixed"] - 3 / 16) < 1e-6 and abs(c["chi_inv"]) < 1e-6


def test_twobody():
    for M in (0.1, 1.0, 1836.0):
        d = 1e-4
        F = T.coulomb_fidelity(M * np.exp(-d), M * np.exp(d))
        assert abs(2 * (1 - F) / (2 * d) ** 2 - 0.75 / (M + 1) ** 2) < 1e-6
        F = T.harmonic_fidelity(M * np.exp(-d), M * np.exp(d))
        assert abs(2 * (1 - F) / (2 * d) ** 2 - T.harmonic_chi_fixed(M)) < 1e-6
        # dilation by s = mu1/mu2 maps the Coulomb state exactly: numerical quadrature
        M2 = 3 * M
        b1, b2 = T.mu(M), T.mu(M2)
        s = b1 / b2
        num = quad(lambda r: 4 * np.pi * r ** 2 * np.exp(-b1 * r) * s ** 1.5 * np.exp(-b2 * s * r), 0, 200 / b1)[0]
        assert abs(num * (b1 ** 3 * b2 ** 3) ** .5 / np.pi - 1) < 1e-9


# ------------------------------------------------------------ entanglement
@pytest.mark.parametrize("M", MASSES)
def test_gaussian_bipartitions_vs_grid_svd(M):
    p = A.params(M)
    a, b = p["a"], p["b"]
    Bq = b * (2 * M + 1) ** 2 / 4
    mats = {"R|s": (a + b / 4, -b / 2, b), "x1|x3": (4 * a, 2 * a / M, (a + Bq) / M ** 2)}
    if M < 5:
        mats["x1|x2"] = (a + Bq, Bq - a, a + Bq)
    for w, Am in mats.items():
        P1, S1 = N.gaussian_1d_purity_grid(*Am, n=600)
        e = A.entanglement(M, which=w)
        assert abs(e["purity"] - P1 ** 3) < 1e-9, w
        assert abs(e["S_vN"] - 3 * S1) < 1e-8, w
    assert A.entanglement(M, which="R|rho")["S_vN"] == 0


@pytest.mark.parametrize("M", MASSES)
def test_partial_wave_entanglement(M):
    e0 = B.entanglement_R_s(M, 0.0, n_r=80)
    ea = A.entanglement(M, which="R|s")
    assert abs(e0["S_vN"] - ea["S_vN"]) < 1e-9 and abs(e0["S_L"] - ea["S_L"]) < 1e-11
    e1 = B.entanglement_R_s(M, 1.0, n_r=80)
    assert abs(e1["trace"] - 1) < 1e-8
    en = N.expectations_numeric(M, 1.0)
    r = en["r"]
    f = en["u"] / r / np.sqrt(4 * np.pi)
    pur = N.purity_R_s_quadrature(f, r, r[1] - r[0], B.params(M, 1.0)["b"] / 16)
    assert abs(pur - (1 - e1["S_L"])) < 5e-5


def test_R_s_large_M_limit():
    lim = B.entanglement_R_s_large_M(1.0)
    assert abs(lim["trace"] - 1) < 1e-12
    assert abs(B.entanglement_R_s(1e6, 1.0, n_r=200)["S_vN"] - lim["S_vN"]) < 2e-3


def test_cm_shell_radius_and_dimensionless_curvature():
    M = 10.0
    p = B.params(M, 1.0)
    y = np.linspace(1e-6, 12 / np.sqrt(p["a"]), 6000)
    nm = B.heavy_density_mid(y, M)
    s2 = 1 / (2 * p["b"] * (2 * M + 1) ** 2)
    xs = np.linspace(0.5, 0.8, 301)
    n = N.cm_density_numeric(xs, y, nm, s2)
    assert abs(xs[np.argmax(n)] - B.cm_shell_radius(M, 1.0)) < 2e-3
    assert abs(B.cm_curvature_dimensionless(M, 1.0) - 2 * s2 * B.cm_curvature_coefficient(M, 1.0)) < 1e-12
    # large-M heavy density about CM: shell -> R_e/2
    assert abs(B.cm_shell_radius(1e6, 1.0) - B.bond_length(1.0) / 2) < 2e-3
