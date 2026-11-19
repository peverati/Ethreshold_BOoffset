import os
import tempfile

import numpy as np
import json
from scipy.integrate import simpson

from threebody import PerimetricBasis, solve, helium_like, overlap_between, symmetric_family
from threebody.convergence import extrapolate
from threebody.scan import scan_masses, write_csv

HE_ALPHAS = (2.6562541, 2.6562541, 4.52188892)


def _he(N=12, alphas=HE_ALPHAS):
    return solve(helium_like(2.0), PerimetricBasis(N, alphas, "S"), 2)


def test_known_helium_expectation_values():
    o = _he(16).observables()
    # well-known values (Drake): <r1> 0.929472295, <r12> 1.422070255, <1/r12> 0.945818448,
    # <1/r1^2> 6.017408866, <r12^2> 2.516439313
    assert abs(o["r02"] - 0.929472295) < 1e-7
    assert abs(o["r01"] - 1.422070255) < 1e-7
    assert abs(o["1/r01"] - 0.945818448) < 1e-7
    assert abs(o["1/r02^2"] - 6.017408866) < 2e-6
    assert abs(o["virial_ratio"] - 2.0) < 1e-6
    assert abs(o["E"] - (o["T"] + o["V"])) < 1e-12


def test_energy_components_reproduce_eigenvalue():
    r = _he(10)
    assert abs(r.energy_components()["E"] - r.energies[0]) < 1e-11


def test_pair_distribution_normalization_and_moments():
    r = _he(10)
    for (i, j) in [(0, 1), (0, 2)]:
        assert abs(r.radial_moment_quadrature(i, j, 0) - 1.0) < 1e-12
        assert abs(r.radial_moment_quadrature(i, j, 1) - r.expect_power(i, j, 1)) < 1e-12
        assert abs(r.radial_moment_quadrature(i, j, -1) - r.expect_power(i, j, -1)) < 1e-12
        R = np.linspace(1e-9, 20, 4001)
        g = r.g_of_r(i, j, R)
        assert abs(simpson(4 * np.pi * R**2 * g, x=R) - 1.0) < 1e-8


def test_overlaps_between_bases():
    r1 = _he(12)
    r2 = _he(10, (2.4, 2.4, 4.2))
    assert abs(overlap_between(r1, r1) - 1.0) < 1e-13
    assert abs(overlap_between(r1, r1, 0, 1)) < 1e-10  # orthogonal eigenstates
    ov = overlap_between(r1, r2)
    assert 0.99999 < ov <= 1.0 + 1e-12
    # mixed-scale overlap with equal scales reproduces the overlap matrix exactly
    assert abs(overlap_between(r1, _he(12), normalized=False) - r1.norm()) < 1e-12


def test_save_roundtrip():
    r = _he(8)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "x.npz")
        r.save(p)
        z = np.load(p)
        meta = json.loads(str(z["meta"]))
        assert meta["basis"]["N"] == 8
        assert np.allclose(z["energies"], r.energies)
        assert np.allclose(z["tensors"][0], r.tensor(0))


def test_extrapolate_geometric():
    E = -1 + 0.3 * 0.5 ** np.arange(6)
    Ee, err, rho, flag = extrapolate(E)
    assert abs(Ee + 1) < 1e-14 and abs(rho - 0.5) < 1e-12 and flag == ""


def test_scan_smoke_and_smoothness():
    rows = scan_masses([1.0, 1.01, 1.02], [8, 10, 12], nstates=2, N_opt=8, verbose=False)
    with tempfile.TemporaryDirectory() as d:
        write_csv(os.path.join(d, "s.csv"), rows)
    E = np.array([r["E0"] for r in rows])
    # second difference small compared to first difference (smooth in M)
    assert abs(E[0] - 2 * E[1] + E[2]) < 0.05 * abs(E[1] - E[0])
    assert all(r["n_below_threshold"] == 1 for r in rows)


def test_hellmann_feynman_mass_derivative():
    al, N, M, h = (1.3, 1.3, 1.2), 10, 4.0, 1e-4
    b = PerimetricBasis(N, al, "S")
    from threebody import MatrixSet
    ms = MatrixSet(b)
    e = [solve(symmetric_family(m), b, 1, matrices=ms).energies[0] for m in (M - h, M + h)]
    r = solve(symmetric_family(M), b, 1, matrices=ms)
    hf = -(r.dE_dinverse_mass(0) + r.dE_dinverse_mass(1)) / M**2
    assert abs((e[1] - e[0]) / (2 * h) - hf) < 1e-8
