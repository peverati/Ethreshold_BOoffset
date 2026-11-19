"""Check of the spurious-root guard in run_exact.py on the known J0v1 case (written 2026-11-16).

The original J0v1 stage 3 gave M* = 9.1976224 (tol 1e-12, margin 0) and 9.2172467 (margin 1e-11) for
"compact N=16 + diffuse alpha_s=0.0003", far below every other crossing (~11.5738): a spurious root.
Required: that basis at tol 1e-12 near M = 9.2 is flagged; the same basis at tol 1e-10, N = 14 / 18 with the
same diffuse set, and N = 16 with alpha_s = 1e-4 / 1e-3, all at their recorded crossings, are not.
Finding (first version of this check): the 9.2 rows are NOT a low spurious E_0 (guards (a), (b) do not fire,
E_0(compact) - E_0 ~ 1e-15) but a near-threshold level of the alpha_s=3e-4 basis at tol 1e-12 whose E - E_th is
flat in M and changes between builds (+6e-12 vs -1e-8 at M = 9.1976); guard (c) (crossing slope) catches it.
The alpha_s=3e-4 tol-1e-12 rows at 11.5738 may or may not be flagged (same artefact, build dependent): reported,
not asserted (None).

Run from the project root (needs J0v1_scales.json; ~5 min with 10 threads):
    THREEBODY_BLAS_THREADS=10 python paper1/unequal_masses/check_guard.py
"""
import gc
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_exact as R  # noqa: E402

T = "J0v1"
sc = R.stage_scales(T)
# (N, alpha_s, tol, M, must_flag): masses are the recorded crossings of the original run (J0v1_crossings.csv)
CASES = [
    (16, 3e-4, 1e-12, 9.197622375347665, True),     # spurious margin-0 crossing
    (16, 3e-4, 1e-12, 9.21724669026899, True),      # spurious margin-1e-11 crossing
    (16, 3e-4, 1e-10, 11.573789329451385, False),   # same basis, tol 1e-10
    (16, 1e-4, 1e-12, 11.57375981612933, False),    # genuine best ladder member
    (16, 1e-3, 1e-12, 11.573888704918156, False),   # ladder member alpha_s = 1e-3
    (14, 3e-4, 1e-12, 11.573788750535398, None),    # N - 2 (same fragile diffuse set: not asserted)
    (18, 3e-4, 1e-12, 11.57378932266815, None),     # N + 2 (same fragile diffuse set: not asserted)
]

fails = 0
bases = {}
for N, a, tol, M, must in CASES:
    key = (N, a)
    if key not in bases:
        bases.clear()
        gc.collect()
        bases[key] = R.LeanBasis(0, R.blocks_for(T, sc, N, alpha_s=a), 2.0, tols=(1e-12, 1e-10))
    b = bases[key]
    e = b.energies(M, 2, tol)
    ec = R.compact_basis(T, N).energies(M, 1)[0]
    why = R.guard(T, b, M, tol, N, margin=0.0)
    ok = must is None or bool(why) == must
    fails += not ok
    print(f"{'OK  ' if ok else 'FAIL'} N={N} alpha_s={a:g} tol={tol:g} M={M:.8f}: kept={b.proj[tol]['kept']} "
          f"lam_min={b.proj[tol]['lam_min']:.1e} E0-E_th={e[0] - R.threshold(M, 2.0):+.6e} "
          f"E1-E_th={e[1] - R.threshold(M, 2.0):+.3e} E0c-E0={ec - e[0]:+.3e} -> "
          f"{'flagged: ' + why if why else 'not flagged'} (expected {'-' if must is None else 'flag' if must else 'no flag'})", flush=True)
print("ALL CASES PASS" if not fails else f"{fails} CASE(S) FAILED")
sys.exit(1 if fails else 0)
