"""Short tests (~15 s).  From the project root:
    python3 -m pytest -p no:cacheprovider paper_data/bo_offset/massfunc/test_massfunc.py -q
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pytest

import massfunc as mf

R_E = 1.99719332


def test_large_R_limit_one_half():
    """THEORY.md Sec. 5: m_v, m_r -> 1/2 (separated atoms), approached as ~R^-3."""
    a20, a40 = mf.mass_functions(20.0), mf.mass_functions(40.0)
    for k in ("m_v", "m_r"):
        assert abs(a40[k] - 0.5) < 1e-3
        assert abs(a40[k] - 0.5) < abs(a20[k] - 0.5)
        ratio = (a20[k] - 0.5) / (a40[k] - 0.5)          # R^-3 tail -> ratio 8
        assert 6.0 < ratio < 10.0


def test_basis_convergence_R2():
    vals = [mf.mass_functions(2.0, level=lv) for lv in (-1, 0, 1)]
    for k in ("m_v", "m_r"):
        spread = max(v[k] for v in vals) - min(v[k] for v in vals)
        assert spread < 1e-9, (k, spread)
    alt = mf.mass_functions(2.0, level=0, alpha_scale=0.7)  # different exponent and quadrature nodes
    assert abs(alt["m_v"] - vals[1]["m_v"]) < 1e-9 and abs(alt["m_r"] - vals[1]["m_r"]) < 1e-9
    assert abs(vals[1]["E0_basis"] - vals[1]["E0_exact"]) < 1e-10


def test_rhs_norms_match_adiabatic_D():
    """<f|f> of the right-hand sides, projected on the basis, vs bo_h2plus D(R) pieces."""
    r = mf.mass_functions(2.0)
    dr, dp = mf.reference_norms(2.0)
    assert r["fpar_norm_in_basis"] / dr == pytest.approx(1.0, abs=1e-3)
    assert r["fperp_norm_in_basis"] / dp == pytest.approx(1.0, abs=1e-3)
    assert abs(r["fpar_overlap_phi0"]) < 1e-8


def test_values_at_Re_vs_data():
    """No tuning: m_v(R_e) vs the kappa^6 value 0.3465 +- 0.0003; m_r(R_e) vs J-spacing value 0.0754 +- 0.001."""
    r = mf.mass_functions(R_E)
    assert abs(r["m_v"] - 0.3465) < 0.0003, r["m_v"]
    assert abs(r["m_r"] - 0.0754) < 0.001, r["m_r"]
