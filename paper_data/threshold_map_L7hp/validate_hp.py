"""validate_L with the double-double assembly (HPCoupledBasis).  Written 2026-11-17.

Same test as paper_data/threshold_map/validate_L.py (H2+ M = 1836.15267247, J = L, single block alphas AL, weights
(1,1,2); delta = E_exact - E_BO+ad compared with the quadratic extrapolation of the J <= 4 deltas from the
reference energies; tolerance 1e-7 + 3 |E(N=20) - E(N=16)|), plus: the N sequence 12..20 must be monotone
(variational) at the chosen projection tol.  Usage (from this folder):
    HP_NPROC=20 python validate_hp.py L tol [L tol ...]        e.g.  validate_hp.py 7 1e-14 7 1e-15 6 1e-14
Writes validate_hp.txt (appends)."""
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh

import tb_hp_local  # noqa: F401
from hp_assembly import HPCoupledBasis

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "solvers"))
from threebody import symmetric_family  # noqa: E402
from threebody.coupled import block  # noqa: E402

M_P06 = 1836.15267247
AL = (9.95811007, 9.95811007, 1.42227562)
REF = {0: -0.597139063079392297758, 1: -0.5968737387847130778, 2: -0.5963452054891147, 3: -0.5955576389803092,
       4: -0.59451716906062}
NS = (12, 14, 16, 18, 20)


def e_boad(J):
    from bo_h2plus import solve_levels
    return float(solve_levels(M_P06, J, adiabatic=True, ppw=4.0).levels[0].E)


def main():
    args = sys.argv[1:]
    jobs = [(int(args[i]), float(args[i + 1])) for i in range(0, len(args), 2)]
    out = open(os.path.join(HERE, "validate_hp.txt"), "a")
    P = lambda t: (print(t, flush=True), out.write(t + "\n"), out.flush())
    P(f"# validate_hp.py {time.strftime('%Y-%m-%d %H:%M')}: HP (double-double) assembly; jobs {jobs}; "
      f"HP_NPROC={os.environ.get('HP_NPROC', '20')} HP_TWO_PASS={os.environ.get('HP_TWO_PASS', '0')}")
    d = {J: REF[J] - e_boad(J) for J in range(5)}
    x = np.array([J * (J + 1) for J in range(5)], float)
    c = np.polyfit(x, [d[J] for J in range(5)], 2)
    s = symmetric_family(M_P06)
    ok_all = True
    bases = {}
    for L, tol in jobs:
        bo = e_boad(L)
        pr = float(np.polyval(c, L * (L + 1)))
        E = {}
        for N in NS:
            if (L, N) not in bases:
                t0 = time.time()
                bases.clear()
                bases[(L, N)] = HPCoupledBasis(L, block(N, AL, (1, 1, 2)), "A" if L % 2 else "S")
                bt = time.time() - t0
            else:
                bt = 0.0
            b = bases[(L, N)]
            b._proj = None
            Hm = b.hamiltonian_projected(s, tol)
            E[N] = float(eigh(Hm, eigvals_only=True, subset_by_index=[0, 0])[0])
            P(f"L={L} tol={tol:g} N={N} size={b.size} kept={Hm.shape[0]} lam_min/lam_max={b._proj['lam_min'] / b._proj['lam_max']:+.1e} "
              f"E={E[N]:.12f} E-E_BO+ad={E[N] - bo:+.4e} (build {bt:.0f}s)")
        mono = all(E[NS[i + 1]] <= E[NS[i]] + 1e-12 for i in range(len(NS) - 1))
        dl = E[20] - bo
        tl = 1e-7 + 3 * abs(E[20] - E[16])
        ok = abs(dl - pr) < tl and mono
        ok_all &= ok
        P(f"L={L} tol={tol:g}: delta(N=20) = {dl:+.5e}, extrapolated {pr:+.5e}, diff {dl - pr:+.1e} (tol {tl:.1e}); "
          f"monotone in N: {mono}  -> {'PASS' if ok else 'FAIL'}")
    P("RESULT: " + ("PASS" if ok_all else "FAIL"))


if __name__ == "__main__":
    main()
