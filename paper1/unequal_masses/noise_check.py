"""Build-to-build noise of E_k - E_th at genuine J0v1 crossings (written 2026-11-16).

Rebuilds "compact N=16 + diffuse alpha_s" (alpha_s = 1e-3, 1e-4) at its recorded crossing mass and prints
gap = E_1 - E_th at M*(1 - 1e-3), M*, M*(1 + 1e-3).  Rounding perturbation = BLAS thread count (results are
bit-identical for the same thread count).  Output of the 2026-11-16 run: noise_check_output.jsonl (threads 2, 5,
10, 5).  Run from the project root:
    for th in 2 5 10 5; do THREEBODY_BLAS_THREADS=$th OPENBLAS_NUM_THREADS=$th python paper1/unequal_masses/noise_check.py; done
"""
import json, os, sys, time
sys.path.insert(0, "paper1/unequal_masses"); import run_exact as R
T = "J0v1"; sc = R.stage_scales(T)
CASES = {1e-3: 11.573888704918156, 1e-4: 11.57375981612933}
out = []
for a, Ms in CASES.items():
    t0 = time.perf_counter()
    b = R.LeanBasis(0, R.blocks_for(T, sc, 16, alpha_s=a), 2.0)
    g = {s: R.gap(T, b, Ms * (1 + s)) for s in (-1e-3, 0.0, 1e-3)}
    out.append(dict(threads=os.environ["THREEBODY_BLAS_THREADS"], alpha_s=a, M=Ms, kept=b.proj[1e-12]["kept"],
                    gap_lo=g[-1e-3], gap=g[0.0], gap_hi=g[1e-3], secs=round(time.perf_counter() - t0)))
    print(json.dumps(out[-1]), flush=True)
    del b
