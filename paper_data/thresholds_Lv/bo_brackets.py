"""BO and BO+ad critical masses that bracket the exact searches of thresholds_Lv (cheap, ~1-3 min).

Written 2026-11-15.  Run from the project root:
    python3 paper_data/thresholds_Lv/bo_brackets.py

1sigma_g sector (paper_data/bo_offset/coupling.py, zero-energy coupling problem, lam = M + 1/2, any L;
c = 0 BO, c = 1 BO+ad counted at its own asymptote):
    (J=1, v=0), (J=1, v=1), (J=3, v=0), (J=4, v=0)   [+ J=2 v=0 as a check against bo_critical_J.txt]
Cross-check of the BO values: bo_h2plus count_levels bisection (Numerov, ppw 16 and 32, Richardson).

2p sigma_u sector (exchange-SYMMETRIC odd L / antisymmetric even L of {X+,X+,e-}): BO only.
U_u(R) = E_el(R, '2pu') + 1/R from bo_h2plus.electronic (unit mass), tabulated and splined;
asymptote -1/2 - 9/(4R^4) (same polarization tail as 1s sigma_g; the g/u splitting is exponentially
small).  Zero-energy shooting gives lam_c for J = 0 and J = 1 (number of bound levels vs lam).
No adiabatic correction for 2p sigma_u (bo_h2plus.adiabatic only implements 1s sigma_g).

Output: bo_brackets.json, bo_brackets.txt.
"""
import json
import math
import os
import sys
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("BO_H2PLUS_THREADS", "4")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "paper_data", "bo_offset"))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.interpolate import CubicSpline  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import coupling as C  # noqa: E402

OUT = []


def P(s=""):
    print(s, flush=True)
    OUT.append(s)


def lam_roots(L, c, lo, hi, n):
    return [float(x) for x in C.critical(L, c, lo, hi, n=n)]


def bo_h2plus_check(J, n, lo, hi):
    from bo_h2plus import count_levels

    def crit(ppw, lo, hi, rtol=1e-11):
        cnt = lambda M: count_levels(M, J, adiabatic=False, ppw=ppw)
        assert cnt(lo) < n <= cnt(hi), (J, n, lo, hi, cnt(lo), cnt(hi))
        while hi / lo - 1 > rtol:
            mid = math.sqrt(lo * hi)
            if cnt(mid) >= n:
                hi = mid
            else:
                lo = mid
        return math.sqrt(lo * hi)
    M1 = crit(16.0, lo, hi)
    M2 = crit(32.0, M1 * 0.99, M1 * 1.01)
    return M2 + (M2 - M1) / 15, abs(M2 - M1) / 15


# ------------------------------------------------------------------ 2p sigma_u
def u_potential():
    from bo_h2plus.electronic import electronic_energy
    R = np.concatenate([np.geomspace(0.05, 2.0, 40, endpoint=False), np.linspace(2.0, 30.0, 141),
                        np.geomspace(30.0, 120.0, 30)[1:]])
    E = np.array([electronic_energy(float(r), state="2pu") for r in R])
    dU = E + 1.0 / R + 0.5          # unit-mass BO potential minus asymptote
    spl = CubicSpline(np.log(R), R * dU)
    RB = R[-1]

    def f(r):
        r = np.asarray(r, float)
        return np.where(r <= RB, spl(np.log(np.clip(r, R[0], RB))) / r, -9 / (4 * r ** 4) - 15 / (2 * r ** 6))
    tail = dU[-1] * RB ** 4
    return R, dU, f, tail


def shoot_u(lam, L, f, RM=4.0e3):
    """Zero-energy solution of -u'' + [L(L+1)/r^2 + lam dU_u(r)] u = 0; returns the growing amplitude
    (J = 0: Wronskian with the zero-energy -C4/r^4 solution; J > 0: coefficient of r^(L+1))."""
    def rhs(r, y):
        return [y[1], (L * (L + 1) / r ** 2 + lam * f(r)) * y[0]]
    r0 = 0.05
    # repulsive core (dU ~ +1/r): start with u ~ r^(L+1)
    y0 = [r0 ** (L + 1), (L + 1) * r0 ** L]
    s = solve_ivp(rhs, (r0, RM), y0, method="DOP853", rtol=1e-11, atol=1e-300)
    u, up = s.y[0, -1], s.y[1, -1]
    nrm = np.max(np.abs(s.y[0]))
    r = RM
    if L == 0:
        b = math.sqrt(lam * 9 / 4)
        th = math.atan2(u / r, (u - r * up) / b)
        amp = math.hypot(u / r, (u - r * up) / b)
        return amp * math.sin(th - b / r) / nrm
    return (L * u + r * up) / ((2 * L + 1) * r ** (L + 1)) / nrm


def u_criticals(f, L, lams):
    F = [shoot_u(x, L, f) for x in lams]
    out = []
    for i in range(len(lams) - 1):
        if np.sign(F[i]) != np.sign(F[i + 1]):
            out.append(float(brentq(shoot_u, lams[i], lams[i + 1], args=(L, f), xtol=1e-9, rtol=1e-13)))
    return out


def main():
    t0 = time.time()
    res = {}
    P(f"# BO / BO+ad critical masses for thresholds_Lv  ({time.strftime('%Y-%m-%d %H:%M')})")
    P("# M_c = lam_c - 1/2 (electron masses); 1s sigma_g sector from paper_data/bo_offset/coupling.py")
    jobs = {1: (1.0, 40.0, 80), 2: (5.0, 20.0, 40), 3: (10.0, 45.0, 60), 4: (20.0, 70.0, 60)}
    for L, (lo, hi, n) in jobs.items():
        for c, name in ((0, "BO"), (1, "BO+ad")):
            lams = lam_roots(L, c, lo, hi, n)
            res[f"g L={L} {name}"] = [x - 0.5 for x in lams]
            P(f"1sg L={L} {name:6s} M_c(v=0,1,..) = " + ", ".join(f"{x - 0.5:.8f}" for x in lams))
    P("# cross-check (pure BO) with bo_h2plus count_levels bisection, ppw 16/32 Richardson")
    for J, n in ((1, 2), (3, 1), (4, 1)):
        Mc = res[f"g L={J} BO"][n - 1]
        v, e = bo_h2plus_check(J, n, Mc * 0.97, Mc * 1.03)
        res[f"g L={J} v={n - 1} BO bo_h2plus"] = v
        P(f"  J={J} v={n - 1}: bo_h2plus {v:.8f} (+/- {e:.1e});  coupling {Mc:.8f};  diff {v - Mc:+.1e}")
    # 2p sigma_u
    P("# 2p sigma_u sector, BO only (no adiabatic correction available)")
    R, dU, f, tail = u_potential()
    i = int(np.argmin(dU))
    P(f"  2pu well (unit mass): min dU = {dU[i]:.6e} Eh at R = {R[i]:.2f} bohr (tabulated grid); "
      f"dU*R^4 at R={R[-1]:.0f}: {tail:.4f} (polarization -9/4 = -2.25)")
    lams = np.geomspace(50.0, 2.0e5, 120)
    for L in (0, 1, 2):
        cr = u_criticals(f, L, lams)
        res[f"u L={L} BO"] = [x - 0.5 for x in cr]
        P(f"2pu L={L} BO     M_c(v=0,1,..) = " + (", ".join(f"{x - 0.5:.6f}" for x in cr[:4]) or "none below M=2e5"))
    res["note_2pu"] = ("2p sigma_u J=0 <-> exchange-antisymmetric L=0; J=1 <-> exchange-SYMMETRIC L=1 ('S'); "
                       "1s sigma_g J odd <-> 'A', J even <-> 'S'.")
    json.dump(res, open(os.path.join(HERE, "bo_brackets.json"), "w"), indent=1)
    P(f"# done in {time.time() - t0:.0f}s")
    open(os.path.join(HERE, "bo_brackets.txt"), "w").write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
