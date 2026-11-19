"""Consistency check of the general-L angular code for L >= 5 (not validated before) on H2+ rotational levels.

Written 2026-11-16 (threshold_map).  The tensor code of solvers/threebody (thresholds_Lv extension) was validated for
L = 3 (Tang et al. 2014, 1e-11) and L = 4 (few 1e-9).  For L = 5, 6, 7 no reference is in hand.  Check used here:
the non-adiabatic remainder delta(J) = E_exact(v=0, J) - E_BO+ad(v=0, J) (BO+ad from solvers/bo_h2plus, DVR) is a
smooth function of x = J(J+1).  delta_ref(J) for J = 0..3 uses the reference energies of Tang et al. 2014 (as in
solvers/threebody/examples/rotational_validation_L34.py) and for J = 4 the converged value of the validated L = 4 code
(thresholds_Lv/validation_h2p_L34.txt, extrapolated -0.59451716906062); it is almost exactly linear in x
(-5.3e-9 per unit x).  A quadratic fit in x is extrapolated to J = 5..7 and compared with delta from the new term lists
(N = 16 and 20; the N-difference estimates the basis error, ~3e-8 at J = 3, 4 for N = 20).  An error in the new term
lists would change E(J) by far more than TOL_DELTA (the rotational energy at J = 5..7 is ~4e-3 to 7e-3 Eh).
(The S+A = unsymmetrized check of thresholds_Lv is NOT repeated: for L >= 5 the unsymmetrized basis needs term lists
for every component pair, another ~1 h of sympy per L; the code path is the one tested for L = 3, 4.)

Usage (project root):  python3 paper_data/threshold_map/validate_L.py 5 [6 7]
Writes validate_L<L>.txt for a single L (validate_L.txt for several; appends); exit code 1 if a check fails (tolerance TOL_DELTA, TOL_UNION).
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
os.environ.setdefault("THREEBODY_BLAS_THREADS", "8")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "solvers"))
import tb_local  # noqa: E402,F401
import numpy as np  # noqa: E402

from threebody import symmetric_family  # noqa: E402
from threebody.coupled import CoupledBasis, block  # noqa: E402

M_P06 = 1836.15267247
AL = (9.95811007, 9.95811007, 1.42227562)     # thresholds_Lv validation scales (L = 1 optimum)
N = 20
TOL_DELTA = 1e-7      # |delta(J) - extrapolated| (Eh); the J <= 4 fit itself is checked below
TOL_UNION = 1e-9


def e_exact(J, n=N):
    b = CoupledBasis(J, block(n, AL, (1, 1, 2)), "A" if J % 2 else "S")
    return float(b.solve(symmetric_family(M_P06), 1, 1e-12, vectors=False)[0][0]), b.size


def e_boad(J):
    from bo_h2plus import solve_levels
    lv = solve_levels(M_P06, J, adiabatic=True, ppw=4.0)
    return float(lv.levels[0].E)


def union_check(L):
    s = symmetric_family(M_P06)
    bl = block(6, (1.2, 1.2, 1.0))
    eS = CoupledBasis(L, bl, "S").solve(s, 4, 1e-12, vectors=False)[0]
    eA = CoupledBasis(L, bl, "A").solve(s, 4, 1e-12, vectors=False)[0]
    eN = CoupledBasis(L, bl, None).solve(s, 8, 1e-12, vectors=False)[0]
    return float(np.max(np.abs(np.sort(np.concatenate([eS, eA]))[:6] - np.sort(eN)[:6])))


def main():
    Ls = [int(x) for x in sys.argv[1:]] or [5]
    out = open(os.path.join(HERE, "validate_L.txt" if len(Ls) > 1 else f"validate_L{Ls[0]}.txt"), "a")
    P = lambda t: (print(t, flush=True), out.write(t + "\n"), out.flush())
    P(f"# validate_L.py {time.strftime('%Y-%m-%d %H:%M')}: H2+ M={M_P06}, alphas {AL}, weights (1,1,2), N={N}; L = {Ls}")
    ok = True
    REF = {0: -0.597139063079392297758, 1: -0.5968737387847130778, 2: -0.5963452054891147, 3: -0.5955576389803092,
           4: -0.59451716906062}
    d = {}
    for J in range(5):
        d[J] = REF[J] - e_boad(J)
        P(f"J={J}: E_ref = {REF[J]:.13f}  E_BO+ad = {REF[J] - d[J]:.13f}  delta_ref = {d[J]:+.5e}")
    x = np.array([J * (J + 1) for J in range(5)], float)
    c = np.polyfit(x, [d[J] for J in range(5)], 2)
    fit_res = max(abs(np.polyval(c, x[i]) - d[i]) for i in range(5))
    P(f"quadratic fit of delta_ref(x), J=0..4: coefficients {c}, max residual {fit_res:.1e}")
    for L in Ls:
        t0 = time.time()
        ex16, n16 = e_exact(L, 16)
        ex, n = e_exact(L)
        bo = e_boad(L)
        dl = ex - bo
        pr = float(np.polyval(c, L * (L + 1)))
        tol = TOL_DELTA + 3 * abs(ex - ex16)
        P(f"L={L}: E_exact N=16 (size {n16}) {ex16:.12f}, N={N} (size {n}) {ex:.12f}; E_BO+ad = {bo:.12f}; "
          f"delta = {dl:+.5e}  extrapolated {pr:+.5e}  diff {dl - pr:+.1e} (tol {tol:.1e})  "
          f"{'OK' if abs(dl - pr) < tol else 'FAIL'}  ({time.time() - t0:.0f}s)")
        ok &= abs(dl - pr) < tol
    P("RESULT: " + ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
