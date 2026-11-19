"""Exact (perimetric-Laguerre, solvers/threebody) energies for the King-Longford-Cox 2013 systems.

Usage (from this folder, python3, THREEBODY_BLAS_THREADS=3):
    python compute.py table     # every distinct mass ratio of their Table III (restartable)
    python compute.py fine      # dense M grid 0.30-0.80 for the g(M) minimum (restartable)
    python compute.py minimum   # root of dg/dM = 0 from exact Hellmann-Feynman derivatives
    python compute.py extra     # the 7 extra a3 points of their fit (a3 = 0.17..0.23)
    python compute.py lit       # Table II systems at the masses of their 'literature' column
then python analyze.py.

Each stage appends rows to its CSV and skips what is already there.
"""
import csv
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq

from common import HERE, read_table3, Eth_reduced, g_of

from threebody import (ThreeBodySystem, symmetric_family, PerimetricBasis, MatrixSet, solve,
                       optimize_alphas, extrapolate)


def log(msg):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    with open(os.path.join(HERE, "run_log.txt"), "a") as fh:
        fh.write(line + "\n")


class Appender:
    def __init__(self, path, fields):
        new = not os.path.exists(path)
        self.fh = open(path, "a", newline="")
        self.w = csv.DictWriter(self.fh, fieldnames=fields, extrasaction="ignore")
        if new:
            self.w.writeheader()

    def row(self, d):
        self.w.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in d.items()})
        self.fh.flush()


def done(path, key):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {r[key] for r in csv.DictReader(fh)}


# ============================================================================ table
def settings(M):
    """Basis settings by mass ratio (as in paper1/mass_scans, one extra N for M < 150)."""
    if M < 150:
        return dict(weights=(1, 1, 1), N_opt=12, Ns=[16, 20, 24, 28, 32])
    if M < 3000:
        return dict(weights=(1, 1, 2), N_opt=20, Ns=[28, 32, 36, 40])
    if M < 5000:
        return dict(weights=(1, 1, 2), N_opt=20, Ns=[32, 36, 40, 44])
    return dict(weights=(1, 1, 2), N_opt=20, Ns=[36, 40, 44, 48])


def system_for(M):
    """Dimensionless system with the smaller of (m1, m3) set to 1.

    M >= 1: masses (M, M, 1);  0 < M < 1: (1, 1, 1/M);  M = 0: (1, 1, inf).
    Returns (system, m_ref_is) where m_ref_is = 'm3' or 'm1' says which physical mass is the unit.
    """
    if M >= 1:
        return symmetric_family(M), "m3"
    m3 = math.inf if M == 0 else 1.0 / M
    return ThreeBodySystem((1.0, 1.0, m3), (1.0, 1.0, -1.0), identical=(0, 1), name=f"11X M={M:g}"), "m1"


# masses used by the "literature" column of King et al. Table II (their Table I caption), units of m_e
LIT = {
    "H-@lit": ("e", 1836.152701),            # Ref. 33: m_p = 1836.152701
    "Mu-@lit": ("e", 206.768262),            # Ref. 40: m_mu = 206.768262
    "muH-@lit": (206.768262, 1836.1526675),  # Ref. 39
    "muD2+@lit": (3670.4829652, 206.768275), # Ref. 41 (Frolov: m_mu = 206.768275, m_d as in this work)
    "Mu2+@lit": (206.768262, "e"),           # Ref. 30
    "H2+@lit": (1836.152701, "e"),           # Ref. 42
    "D2+@lit": (3670.483014, "e"),           # Ref. 43
}


def lit_items():
    out = []
    for k, (m1, m3) in LIT.items():
        m1 = 1.0 if m1 == "e" else m1
        m3 = 1.0 if m3 == "e" else m3
        out.append((k, m1 / m3))
    return out


def extra_items():
    """King et al. added 7 points between (tau,tau,t) (a3 = 0.2403) and (p,p,t) (a3 = 0.1431) in steps of
    0.01 in a3 for their fit; the exact a3 values are not printed. We take a3 = 0.17 ... 0.23."""
    return [(f"a3={a:.2f}", 2 * a / (1 - a)) for a in np.round(np.arange(0.17, 0.2301, 0.01), 2)]


def stage_table(items=None, fname="ours_table.csv"):
    path = os.path.join(HERE, fname)
    fields = ["key", "M", "unit", "weights", "N_opt", "alpha_12", "alpha_3", "Ns", "sizes", "E_seq",
              "E_var", "E_ext", "E_err", "E_rho", "E_flag", "E_th", "virial", "time_s"]
    out = Appender(path, fields)
    have = done(path, "key")
    if items is None:
        items = [(f"{M:.12g}", M) for M in sorted({round(r["M"], 12) for r in read_table3()})]
    for key, M in items:
        if key in have:
            continue
        t0 = time.time()
        st = settings(M)
        s, unit = system_for(M)
        guess = None
        if M > 150:
            guess = (8.0 * (M / 1836.15) ** 0.25,) * 2 + (1.4,)
        al, _, nev = optimize_alphas(s, st["N_opt"], "S", x0=guess, weights=st["weights"])
        seq, sizes, last = [], [], None
        for N in st["Ns"]:
            b = PerimetricBasis(N, al, "S", st["weights"])
            r = solve(s, b, nstates=1, matrices=MatrixSet(b))
            seq.append(float(r.energies[0]))
            sizes.append(b.size)
            last = r
            log(f"table M={key} N={N} size={b.size} E={seq[-1]:.15f} dropped={r.info['dropped']}")
        ext = extrapolate(seq)
        comp = last.energy_components(0)
        out.row(dict(key=key, M=M, unit=unit, weights=json.dumps(st["weights"]), N_opt=st["N_opt"],
                     alpha_12=al[0], alpha_3=al[2], Ns=json.dumps(st["Ns"]), sizes=json.dumps(sizes),
                     E_seq=json.dumps(seq), E_var=seq[-1], E_ext=ext[0], E_err=ext[1], E_rho=ext[2],
                     E_flag=ext[3], E_th=s.threshold(), virial=comp["virial_ratio"], time_s=time.time() - t0))
        log(f"table M={key} done: E_ext={ext[0]:.15f} err={ext[1]:.1e} flag={ext[3]} t={time.time() - t0:.0f}s")


# ============================================================================ fine grid
FINE_LO, FINE_HI, FINE_STEP = 0.30, 0.80, 0.0025
PANEL = 0.05
FINE_NS = (28, 32)


def g_and_dg(E, dE, M):
    Eth = Eth_reduced(M)
    dEth = -1.0 / (2.0 * (M + 1.0) ** 2)
    return E / Eth - 1.0, dE / Eth - E * dEth / Eth ** 2


def point(ms, b, M):
    r = solve(symmetric_family(M), b, nstates=1, matrices=ms)
    E = float(r.energies[0])
    dE = -(r.dE_dinverse_mass(0, 0) + r.dE_dinverse_mass(1, 0)) / M ** 2   # exact at fixed basis
    return E, dE


def stage_fine():
    path = os.path.join(HERE, "fine_grid.csv")
    fields = ["M", "panel_center", "alpha_X", "alpha_e"]
    for N in FINE_NS:
        fields += [f"E_N{N}", f"dEdM_N{N}", f"g_N{N}", f"dgdM_N{N}"]
    out = Appender(path, fields)
    have = {round(float(k), 6) for k in done(path, "M")}
    grid = np.round(np.arange(FINE_LO, FINE_HI + 1e-9, FINE_STEP), 6)
    centers = np.round(np.arange(FINE_LO, FINE_HI + 1e-9, PANEL), 6)
    for c in centers:
        pts = [M for M in grid if abs(M - c) <= PANEL / 2 + 1e-9 and M not in have
               and (abs(M - c) < PANEL / 2 - 1e-9 or M < c)]   # each boundary point to the upper panel
        if not pts:
            continue
        t0 = time.time()
        al, _, _ = optimize_alphas(symmetric_family(float(c)), 12, "S")
        mats = {}
        for N in FINE_NS:
            b = PerimetricBasis(N, al, "S")
            mats[N] = (MatrixSet(b), b)
        log(f"fine panel c={c:.3f} alphas={al[0]:.6f},{al[2]:.6f} npts={len(pts)} build {time.time() - t0:.0f}s")
        for M in pts:
            row = dict(M=float(M), panel_center=float(c), alpha_X=al[0], alpha_e=al[2])
            for N in FINE_NS:
                E, dE = point(*mats[N], float(M))
                g, dg = g_and_dg(E, dE, float(M))
                row.update({f"E_N{N}": E, f"dEdM_N{N}": dE, f"g_N{N}": g, f"dgdM_N{N}": dg})
            out.row(row)
        log(f"fine panel c={c:.3f} done t={time.time() - t0:.0f}s")


# ============================================================================ minimum
def stage_minimum():
    """Root of dg/dM (exact HF derivative of the variational eigenvalue at fixed basis).

    Alpha sets optimized at M_a (N_opt 12 and 16); N = 20..32 nested; also alphas from 0.40 and 0.60
    to show the result does not depend on the scale choice."""
    path = os.path.join(HERE, "minimum.csv")
    fields = ["M_alpha", "N_opt", "alpha_X", "alpha_e", "N", "size", "M_min", "g_min", "E_at_min",
              "d2g_dM2", "n_eval", "time_s"]
    out = Appender(path, fields)
    have = set()
    if os.path.exists(path):
        with open(path) as fh:
            have = {(r["M_alpha"], r["N_opt"], r["N"]) for r in csv.DictReader(fh)}
    cases = [(0.49, 12), (0.49, 16), (0.40, 12), (0.60, 12)]
    for Ma, Nopt in cases:
        Ns = [20, 24, 28, 32] if (Ma, Nopt) == (0.49, 12) or (Ma, Nopt) == (0.49, 16) else [28, 32]
        al = None
        for N in Ns:
            if (repr(Ma), str(Nopt), str(N)) in have:
                continue
            if al is None:
                al, _, _ = optimize_alphas(symmetric_family(Ma), Nopt, "S")
            t0 = time.time()
            b = PerimetricBasis(N, al, "S")
            ms = MatrixSet(b)
            cnt = [0]

            def dg(M):
                cnt[0] += 1
                E, dE = point(ms, b, M)
                return g_and_dg(E, dE, M)[1]

            Mmin = brentq(dg, 0.40, 0.60, xtol=1e-12, rtol=1e-14)
            E, dE = point(ms, b, Mmin)
            g = g_and_dg(E, dE, Mmin)[0]
            h = 1e-3
            d2 = (dg(Mmin + h) - dg(Mmin - h)) / (2 * h)
            out.row(dict(M_alpha=Ma, N_opt=Nopt, alpha_X=al[0], alpha_e=al[2], N=N, size=b.size, M_min=Mmin,
                         g_min=g, E_at_min=E, d2g_dM2=d2, n_eval=cnt[0], time_s=time.time() - t0))
            log(f"minimum alphas@{Ma} Nopt={Nopt} N={N}: M_min={Mmin:.12f} g_min={g:.12f} t={time.time() - t0:.0f}s")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "table"
    {"table": stage_table, "fine": stage_fine, "minimum": stage_minimum,
     "lit": lambda: stage_table(lit_items(), "ours_lit.csv"),
     "extra": lambda: stage_table(extra_items(), "ours_extra.csv")}[stage]()
