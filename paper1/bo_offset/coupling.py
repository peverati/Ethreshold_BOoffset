"""Critical masses as coupling-constant eigenvalues (zero-energy problem).

BO (+ad) radial equation for {X+,X+,e-}, scaled to unit-mass electronic units r = mu_e R:
    -u'' + L(L+1)/r^2 u + lam*dU(r) u + c*dD(r) u = 0,   lam = M/mu_e = M + 1/2 (exact identity)
dU = U(r) + 1/2 (unit-mass BO potential minus its asymptote), dD = D(r) - 1/4 (adiabatic, minus asymptote).
c = 0: BO,  c = 1: BO + adiabatic (counted at its own asymptote).  M_c = lam_c - 1/2.
Critical lam: zero-energy solution with no growing part (L=0: u -> const; L>0: u ~ r^-L).
Also returns first-order weights used for the mass-correction argument.
Run from repo root:  python3 paper1/bo_offset/coupling.py
"""
import os, sys
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "solvers"))
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from bo_h2plus import potential as P

R0, RB, RM = 1e-4, 200.0, 2.0e4
rg = np.geomspace(R0, RB, 40000)
_rdU = CubicSpline(np.log(rg), rg * (P.U(rg) + 0.5))      # r*dU smooth (1/r at 0)
_dD = CubicSpline(np.log(rg), P.D(rg) - 0.25)
DHI = P.D(np.array([RB]))[0] - 0.25

def dU(r):
    r = np.asarray(r, float)
    return np.where(r <= RB, _rdU(np.log(np.minimum(r, RB))) / r, -9 / (4 * r**4) - 15 / (2 * r**6))
def dD(r):
    r = np.asarray(r, float)
    return np.where(r <= RB, _dD(np.log(np.minimum(r, RB))), DHI * (RB / r)**4)

def shoot(lam, L, c, dense=False):
    def f(r, y):
        return [y[1], (L * (L + 1) / r**2 + lam * dU(r) + c * dD(r)) * y[0]]
    r0 = 1e-3
    a1 = lam * 1.0 / (2 * L + 2)          # Frobenius: dU ~ 1/r at small r
    y0 = [r0**(L + 1) * (1 + a1 * r0), (L + 1) * r0**L + (L + 2) * a1 * r0**(L + 1)]
    s = solve_ivp(f, (r0, RM), y0, method="DOP853", rtol=1e-12, atol=1e-30, dense_output=dense)
    u, up = s.y[0, -1], s.y[1, -1]
    r = RM
    if L == 0:
        C4 = lam * 9 / 4 - c * DHI * RB**4
        b = np.sqrt(C4)
        th = np.arctan2(u / r, (u - r * up) / b)          # u = a r sin th, u - r u' = a b cos th
        amp = np.hypot(u / r, (u - r * up) / b)
        A = amp * np.sin(th - b / r)
    else:
        A = (L * u + r * up) / ((2 * L + 1) * r**(L + 1))
    norm = np.max(np.abs(s.y[0]))
    return (A / norm, s) if dense else A / norm

def critical(L, c, lam_lo, lam_hi, n=60):
    lams = np.linspace(lam_lo, lam_hi, n)
    F = np.array([shoot(l, L, c) for l in lams])
    roots = []
    for i in range(n - 1):
        if np.sign(F[i]) != np.sign(F[i + 1]):
            roots.append(brentq(shoot, lams[i], lams[i + 1], args=(L, c), xtol=1e-13, rtol=1e-14))
    return roots

def weights(lam, L, c):
    """Integrals for first-order perturbation theory at the critical lam."""
    _, s = shoot(lam, L, c, dense=True)
    r = np.concatenate([np.geomspace(1e-3, 50, 200001), np.geomspace(50, RM, 200001)[1:]])
    y = s.sol(r); u, up = y[0], y[1]
    u = u / np.max(np.abs(u)); up = up / np.max(np.abs(y[0]))
    tr = lambda g: np.trapz(g, r)
    B = tr(-dU(r) * u**2)
    I = dict(B=B, Drad=tr(up**2), Dcen=tr(L * (L + 1) * u**2 / r**2), Dad=tr(dD(r) * u**2), lam=lam)
    return I, r, u, up

if __name__ == "__main__":
    out = os.path.join(here, "coupling_results.csv")
    rows = []
    for L, lo, hi in [(0, 0.5, 64.0), (1, 3.0, 8.0), (2, 9.0, 15.0)]:
        for c in (0, 1):
            for v, lam in enumerate(critical(L, c, lo, hi)):
                I, r, u, up = weights(lam, L, c)
                rows.append((L, v, c, lam, lam - 0.5, I["B"], I["Drad"], I["Dcen"], I["Dad"]))
                print(rows[-1], flush=True)
    with open(out, "w") as fh:
        fh.write("# zero-energy coupling eigenproblem; lam = M + 1/2; c=0 BO, c=1 BO+ad\n")
        fh.write("# integrals with u normalized to max|u|=1 (ratios only meaningful): B=int(-dU)u^2, Drad=int u'^2, Dcen=int L(L+1)u^2/r^2, Dad=int dD u^2\n")
        fh.write("L,v_index,c,lam_c,M_c,B,Drad,Dcen,Dad\n")
        for rw in rows:
            fh.write(",".join(str(x) for x in rw) + "\n")
    print("wrote", out)
