"""Mass scans and BO critical masses for {X+, X+, e-}."""
from __future__ import annotations

import csv
import math

import numpy as np

from . import potential as P
from .electronic import equilibrium
from .nuclear import solve_levels, count_levels, thresholds

HARTREE_CM = 219474.6313632  # CODATA 2018 hartree in cm^-1
M_H2PLUS = 1836.15267343     # proton/electron mass ratio (CODATA 2018)

_EQ1 = None


def _eq1():
    global _EQ1
    if _EQ1 is None:
        _EQ1 = equilibrium()
    return _EQ1


def m_eff(M):
    """Effective mass of the scaled BO problem: mu_N/mu_e = (2M+1)/4.

    The BO problem for mass M is exactly mu_e times the unit-mass problem
    -(1/(2 m_eff)) d2/dR'2 + U(R') + mu_e J(J+1)/(2 mu_N R'^2) with R' = mu_e R,
    so for J = 0 the level count n(M) depends on M only through m_eff.
    """
    return (2.0 * M + 1.0) / 4.0


COLUMNS = [
    "M", "mu_N", "mu_e", "m_eff",
    "E_thr_BO", "E_thr_exact", "thr_exact_minus_BO",
    "E0", "binding0", "err_E0", "n_levels_J0", "n_levels_J0_dvr",
    "E1_minus_E0", "E1_minus_E0_cm",
    "EJ1_minus_EJ0", "EJ2_minus_EJ0", "EJ1_minus_EJ0_cm", "EJ2_minus_EJ0_cm",
    "R_mean", "dR", "lindemann", "inv_R2", "B0_from_invR2_cm", "T_kin0",
    "R_e", "E_min", "k", "omega_harm", "omega_harm_cm", "dR_harm", "E0_harm",
    "E0_adiabatic", "E_thr_BO_adiabatic",
]


def scan_row(M, adiabatic=False):
    """One CSV row for mass M (pure BO, plus BO+adiabatic E0 if adiabatic=True)."""
    me = P.mu_e(M)
    muN = M / 2.0
    thr_bo, thr_ex = thresholds(M)
    eq1 = _eq1()
    Re, Emin, k = eq1.Re / me, me * eq1.Emin, me ** 3 * eq1.k
    om = math.sqrt(k / muN)
    row = dict(M=M, mu_N=muN, mu_e=me, m_eff=m_eff(M), E_thr_BO=thr_bo, E_thr_exact=thr_ex,
               thr_exact_minus_BO=thr_ex - thr_bo, R_e=Re, E_min=Emin, k=k,
               omega_harm=om, omega_harm_cm=om * HARTREE_CM,
               dR_harm=(2 * muN * om) ** -0.5, E0_harm=Emin + om / 2)
    nan = float("nan")
    s0 = solve_levels(M, 0)
    row["n_levels_J0"] = count_levels(M, 0)
    row["n_levels_J0_dvr"] = s0.n_dvr
    lv = s0.levels
    if lv:
        l0 = lv[0]
        row.update(E0=l0.E, binding0=l0.binding, err_E0=l0.err_E, R_mean=l0.R_mean, dR=l0.dR,
                   lindemann=l0.lindemann, inv_R2=l0.inv_R2,
                   B0_from_invR2_cm=l0.inv_R2 / (2 * muN) * HARTREE_CM, T_kin0=l0.T_kin)
        d01 = lv[1].E - l0.E if len(lv) > 1 else nan
        row.update(E1_minus_E0=d01, E1_minus_E0_cm=d01 * HARTREE_CM)
        for J in (1, 2):
            sJ = solve_levels(M, J, check=False)
            d = sJ.levels[0].E - l0.E if sJ.levels else nan
            row[f"EJ{J}_minus_EJ0"] = d
            row[f"EJ{J}_minus_EJ0_cm"] = d * HARTREE_CM
    if adiabatic:
        sa = solve_levels(M, 0, adiabatic=True, check=False)
        row["E0_adiabatic"] = sa.levels[0].E if sa.levels else nan
        row["E_thr_BO_adiabatic"] = thresholds(M, adiabatic=True)[0]
    return {c: row.get(c, nan) for c in COLUMNS}


def mass_scan(masses, out=None, adiabatic=False, verbose=False):
    rows = []
    for M in masses:
        r = scan_row(float(M), adiabatic)
        rows.append(r)
        if verbose:
            print(f"M={M:g}: n={r['n_levels_J0']} E0={r['E0']:.12f}")
    if out:
        with open(out, "w", newline="") as fh:
            fh.write("# BO (Born-Oppenheimer) predictions for {X+,X+,e-}; atomic units unless *_cm.\n")
            fh.write("# E_thr_BO = -mu_e/2 (BO limit); E_thr_exact = -M/(2(M+1)) (true X + (Xe) limit); "
                     "they differ by ~1/(4M).\n")
            fh.write("# n_levels_J0: node count at the BO threshold; *_harm: harmonic estimates at R_e of mu_e U(mu_e R).\n")
            w = csv.DictWriter(fh, fieldnames=COLUMNS)
            w.writeheader()
            for r in rows:
                w.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in r.items()})
    return rows


def critical_mass(n, J=0, rtol=1e-11, ppw=8.0, box=3.0, M_lo=1e-3, M_hi=None):
    """BO prediction: smallest M at which the J-rotational problem has n bound levels.

    (i.e. the mass at which level v = n-1 appears at the BO threshold -mu_e/2).
    Bisection in log M on the integer node count.
    """
    cnt = lambda M: count_levels(M, J, ppw=ppw, box=box)
    if cnt(M_lo) >= n:
        raise ValueError("lower bracket already has n levels")
    if M_hi is None:
        M_hi = 1.0
        while cnt(M_hi) < n:
            M_lo = M_hi
            M_hi *= 2.0
            if M_hi > 1e9:
                raise RuntimeError("no bracket below 1e9")
    while M_hi / M_lo - 1 > rtol:
        mid = math.sqrt(M_lo * M_hi)
        if cnt(mid) >= n:
            M_hi = mid
        else:
            M_lo = mid
    return math.sqrt(M_lo * M_hi)


def critical_masses(nmax=6, J=0, ppw=(16.0, 32.0)):
    """BO critical masses M_n (n = 1..nmax): level v = n-1 appears at the BO threshold.

    Bisection on the node count with two Numerov step sizes (dx ratio 2); the
    count converges as O(dx^4), so the Richardson value M_c = M2 + (M2 - M1)/15
    is returned and |M2 - M1|/15 is given as the error estimate of the
    finer, unextrapolated value (the extrapolated one is better than that).
    """
    out = []
    lo = 1e-3
    for n in range(1, nmax + 1):
        M1 = critical_mass(n, J, ppw=ppw[0], M_lo=lo)
        M2 = critical_mass(n, J, ppw=ppw[1], M_lo=M1 * 0.99, M_hi=M1 * 1.01)
        Mc = M2 + (M2 - M1) / 15.0
        out.append(dict(n=n, v_new=n - 1, J=J, M_c=Mc, m_eff=m_eff(Mc), M_c_dx1=M1, M_c_dx2=M2,
                        est_abs_err=abs(M2 - M1) / 15.0))
        lo = Mc
    return out
