"""Shared definitions for the King, Longford & Cox (2013) comparison.

Masses: CODATA 2012, as in Table I of King et al., JCP 139, 224306 (2013), in units of m_e.
All their systems are of the form {m1, m1, m3}; the energy is mass-scaling covariant,
E(m1, m1, m3) = m_ref * E(m1/m_ref, m1/m_ref, m3/m_ref), so each is computed as one
dimensionless {M, M, 1} problem (M = m1/m3) and rescaled.
"""
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "3")

MASS = {  # CODATA 2012 (King et al. Table I), units of m_e
    "e": 1.0,
    "mu": 206.7682843,
    "tau": 3477.15,
    "p": 1836.15267245,
    "d": 3670.4829652,
    "t": 5496.9215267,
    "inf": math.inf,
}


def mass_of(label):
    return MASS[label] if label in MASS else float(label)


def parse_paren(s):
    """'-0.599506(9)' -> (value -0.5995069, last certain decimal place 6, uncertain digits '9')."""
    s = s.strip()
    m = re.fullmatch(r"([+-]?\d+)\.(\d*)\((\d+)\)", s)
    if not m:
        v = float(s)
        dec = len(s.split(".")[1]) if "." in s else 0
        return v, dec, ""
    ip, fp, unc = m.groups()
    return float(f"{ip}.{fp}{unc}"), len(fp), unc


def read_table3():
    with open(os.path.join(HERE, "king2013_table3.csv")) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["m1v"], r["m3v"] = mass_of(r["m1"]), mass_of(r["m3"])
        E, dec, unc = parse_paren(r["E_printed"])
        if r["system"] == "(4.7,4.7,1)" and E > 0:
            E = -E   # missing minus sign in the source
        r["E_their"], r["dec_certain"], r["unc_digits"] = E, dec, unc
        r["M"] = 0.0 if math.isinf(r["m3v"]) else r["m1v"] / r["m3v"]
    return rows


def a3_of(M):
    return M / (M + 2.0)


def M_of_a3(a3):
    return 2.0 * a3 / (1.0 - a3)


def Eth_reduced(M):
    """Threshold of {M, M, 1} in units where the odd particle has mass 1."""
    return -M / (2.0 * (M + 1.0))


def g_of(E_red, M):
    """Relative excess binding g = (E0 - E_th)/E_th, scale invariant."""
    Eth = Eth_reduced(M)
    return (E_red - Eth) / Eth
