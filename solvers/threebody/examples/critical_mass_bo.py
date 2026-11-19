"""BO and BO + adiabatic-correction critical masses for {X+, X+, e-} (J = 0).

Uses only the public API of solvers/bo_h2plus (no change to that package):
  * count_levels(M, 0, adiabatic) : node count at the model's own asymptote
      BO:     -mu_e/2
      BO+ad:  -mu_e/2 + mu_e^2/(4M)   (E_ad(R -> inf) = mu_e^2/(4M))
  * solve_levels(M, 0, adiabatic=True) : level energies (for the variant that
    measures BO+ad levels against the *exact* threshold -M/(2(M+1))).

Critical mass n = mass at which the n-th level appears.  Bisection in M on the
integer count at two Numerov step sizes (ppw 16 and 32); count converges as
O(dx^4) -> Richardson M = M2 + (M2 - M1)/15, error |M2 - M1|/15 (as in bo_h2plus).

Run from solvers/:  python threebody/examples/critical_mass_bo.py OUTDIR
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from scipy.optimize import brentq  # noqa: E402

from bo_h2plus import count_levels, solve_levels  # noqa: E402
from bo_h2plus.scan import critical_masses  # noqa: E402


def crit(n, adiabatic, ppw, lo, hi, rtol=1e-11):
    cnt = lambda M: count_levels(M, 0, adiabatic=adiabatic, ppw=ppw)
    assert cnt(lo) < n <= cnt(hi), (n, lo, hi, cnt(lo), cnt(hi))
    while hi / lo - 1 > rtol:
        mid = math.sqrt(lo * hi)
        if cnt(mid) >= n:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def crit_exact_threshold(n, lo, hi):
    """BO+ad: mass where level v = n-1 reaches the exact threshold -M/(2(M+1))."""
    def f(M):
        lv = solve_levels(M, 0, adiabatic=True, check=False).levels
        Eth = -M / (2 * (M + 1))
        return (lv[n - 1].E - Eth) if len(lv) >= n else 1.0
    return brentq(f, lo, hi, xtol=1e-9)


def main(out):
    lines = ["# BO and BO+adiabatic critical masses, J = 0 (n = number of bound levels after M_c)",
             "# model,n,M_c,est_err,M_c_ppw16,M_c_ppw32"]
    bo = critical_masses(nmax=4)
    for r in bo:
        lines.append(f"BO,{r['n']},{r['M_c']:.10f},{r['est_abs_err']:.1e},{r['M_c_dx1']:.10f},{r['M_c_dx2']:.10f}")
    brackets = {1: (0.3, 3.0), 2: (8.0, 14.0), 3: (25.0, 40.0), 4: (50.0, 75.0)}
    lines.append("# BO+ad: n = 1 is bound already at M = 0.05 (count_levels), so no n = 1 critical mass in the range studied")
    for n in range(2, 5):
        lo, hi = brackets[n]
        M1 = crit(n, True, 16.0, lo, hi)
        M2 = crit(n, True, 32.0, M1 * 0.99, M1 * 1.01)
        Mc = M2 + (M2 - M1) / 15
        lines.append(f"BO+ad(own asymptote),{n},{Mc:.10f},{abs(M2 - M1) / 15:.1e},{M1:.10f},{M2:.10f}")
    for n in range(2, 5):
        lo, hi = brackets[n]
        Mc = crit_exact_threshold(n, lo, hi * 1.3)
        lines.append(f"BO+ad(vs exact threshold),{n},{Mc:.8f},~1e-6 (DVR),,")
    txt = "\n".join(lines) + "\n"
    print(txt)
    with open(os.path.join(out, "bo_critical_masses.csv"), "w") as fh:
        fh.write(txt)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
