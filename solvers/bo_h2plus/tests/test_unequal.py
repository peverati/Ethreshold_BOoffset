"""Tests of bo_h2plus.unequal (added 2026-11-16).  Run from solvers/:  python -m pytest bo_h2plus/tests -q"""
import math
import os

import numpy as np
import pytest

from bo_h2plus import electronic as el
from bo_h2plus import unequal as UQ
from bo_h2plus.adiabatic import D_of_R

HERE = os.path.dirname(os.path.abspath(__file__))
SCHWENKE = os.path.join(HERE, "..", "..", "..", "paper_data", "bo_offset", "massfunc", "literature_data",
                        "schwenke2001_EPAPS_h2p_tab.txt")


@pytest.mark.parametrize("r", [1.0, 1.5, 2.0, 3.0])
def test_lambda_map(r):
    for M in (0.7, 10.0, 1836.0):
        m = UQ.Masses(M, r)
        assert abs(m.lam - UQ.lam_of_M(M, r)) < 1e-12 * m.lam
        assert abs(UQ.M_of_lam(m.lam, r) - M) < 1e-11 * M
    if r == 1.0:
        assert abs(UQ.lam_of_M(10.0, 1.0) - 10.5) < 1e-14


@pytest.mark.parametrize("r", [1.5, 2.0])
def test_thresholds_scaled(r):
    """Exact open/closed thresholds relative to the g BO+ad asymptote -> -beta, +beta (O(1/M) corrections)."""
    b = UQ.beta_of(r)
    for M in (1e4, 1e6):
        lo, hi = UQ.Masses(M, r).thresholds_scaled()
        assert abs(lo + b) < 20.0 / M and abs(hi - b) < 20.0 / M
    m = UQ.Masses(10.0, r)
    assert m.E_th == m.E_th_B < m.E_th_A


def test_fields_R2():
    q = UQ.gu_functions(2.0)
    assert abs(q["norm_g"] - 1) < 1e-10 and abs(q["norm_u"] - 1) < 1e-10 and abs(q["overlap_gu"]) < 1e-12
    assert abs(q["D_g"] - D_of_R(2.0)) < 1e-8
    assert abs(q["t_g"] - UQ.kinetic_t_virial(2.0)) < 1e-7
    assert abs(q["t_u"] - UQ.kinetic_t_virial(2.0, "2pu")) < 1e-7
    # hypervirial <g|d_z|u> = (E_u - E_g) <g|z|u>, antisymmetry, parity zero
    assert abs(q["dz_gu"] - (q["U_u"] - q["U_g"]) * q["z_gu"]) < 1e-8
    assert abs(q["dz_gu"] + q["dz_ug"]) < 1e-9
    assert abs(q["F_overlap_gu"]) < 1e-12
    assert abs(q["U_g"] - el.total_energy(2.0)) < 1e-13


def test_asymptotics_and_localized_eigenvalues():
    q = UQ.gu_functions(30.0)
    assert abs(q["k_gu"] + 1) < 1e-4 and abs(q["D_u"] - 0.25) < 1e-4 and abs(q["t_u"] - 1) < 1e-4
    assert abs(q["dz_gu"]) < 1e-9
    for r in (1.5, 2.0, 3.0):
        b = UQ.beta_of(r)
        K = np.array([[q["D_g"] + b * b * q["t_g"], b * q["k_gu"]], [b * q["k_gu"], q["D_u"] + b * b * q["t_u"]]])
        ev = np.linalg.eigvalsh(K)
        # electron on B: 1/(1+r)^2 ; on A: r^2/(1+r)^2  (atomic reduced-mass shifts)
        assert abs(ev[0] - 1 / (1 + r) ** 2) < 2e-4 and abs(ev[1] - r * r / (1 + r) ** 2) < 2e-4


def _schwenke_hdplus_ad():
    lines = open(SCHWENKE).read().splitlines()
    i = next(k for k, s in enumerate(lines) if s.startswith("Adiabatic correction functions for HD+"))
    out = {}
    for s in lines[i + 2:]:
        p = s.split()
        if len(p) != 3:
            if out:
                break
            continue
        out[float(p[0])] = float(p[1]) * 1e-4
    return out


@pytest.mark.skipif(not os.path.exists(SCHWENKE), reason="literature table not present")
def test_dboc_vs_schwenke_hdplus():
    """Clamped-nuclei DBOC of HD+ = (1/m_p + 1/m_d)(D/2 + t_g/8) (THEORY.md section 3) vs Schwenke 2001 EPAPS
    (aug-cc-pV6Z, uncontracted). His masses are not printed; CODATA-2012-like values are used (m_d quoted from
    memory, irrelevant at this tolerance). The tolerance
    (1e-3 relative) reflects his Gaussian basis, not ours."""
    tab = _schwenke_hdplus_ad()
    mp, md = 1836.15267245, 3670.48296513
    for R in (2.0, 3.0, 5.0):
        q = UQ.gu_functions(R)
        ours = (1 / mp + 1 / md) * (q["D_g"] / 2 + q["t_g"] / 8)
        assert abs(ours / tab[R] - 1) < 1e-3, (R, ours, tab[R])
