"""Rounding-safe upper bounds on M_ck.

The fixed-basis crossing M* (E_k = E_th) is a rigorous upper bound only up to
rounding (~1e-13 hartree in E_k).  Here the crossing is recomputed with a margin,
E_k(M) = E_th(M) - margin (margin = 1e-11 and 1e-10 >> rounding, tol = 1e-12 and
1e-10), and the fixed-basis slope d(E_k - E_th)/dM at M* is reported (Hellmann-
Feynman), which converts the energy rounding into a mass uncertainty.

Run from solvers/threebody:  python examples/critical_mass_bounds.py OUTDIR
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from scipy.optimize import brentq  # noqa: E402

from critical_mass import CONFIG, basis  # noqa: E402
from threebody import family_threshold  # noqa: E402


def main(out):
    lines = ["state,variant,tol,margin,M_bound,slope_at_Mstar"]
    for k in (1, 2, 3):
        c = CONFIG[k]
        for label, kw in [("alpha_s=1e-4", {}), ("alpha_s=1e-4 nd=60 nd2=45", dict(nd=60, nd2=45))]:
            mb = basis(k, 1e-4, **kw)
            for tol in (1e-12, 1e-10):
                f0 = lambda M: mb.energies(M, k + 1, tol)[k] - family_threshold(M)
                Ms = brentq(f0, c["M_lo"], c["M_hi"], xtol=1e-10)
                ob = mb.observables(Ms, state=k, tol=tol)
                slope = ob["dE_dM"] + 1 / (2 * (Ms + 1) ** 2)
                lines.append(f"{k},{label},{tol:g},0,{Ms:.9f},{slope:.3e}")
                for margin in (1e-11, 1e-10):
                    Mm = brentq(lambda M: f0(M) + margin, c["M_lo"], c["M_hi"], xtol=1e-10)
                    lines.append(f"{k},{label},{tol:g},{margin:g},{Mm:.9f},")
                print("\n".join(lines[-3:]), flush=True)
    with open(os.path.join(out, "upper_bounds.csv"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
