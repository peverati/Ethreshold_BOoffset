"""Generate (once) the sympy term lists of natural-parity L sectors into the LOCAL cache (tb_local.py).

Written 2026-11-16 (threshold_map).  One-time cost, single core (sympy): L = 4 took ~15 min in thresholds_Lv;
L = 5, 6, 7 are measured in PROGRESS.md.  Run from the project root:
    python3 paper1/threshold_map/gen_terms.py 6
Building a 1-function basis per L and symmetry triggers every term list the coupled basis needs
(operators S, T01, T2, V01, V02+V12 between all channel pairs c, d and c, Pd).
Each L is generated into its own file _cache/terms_L<L>.pkl (so several L can run in parallel safely) and then
merged into _cache/angular_terms_map.pkl under a file lock.  If the merged cache already holds L, nothing is done.
A per-L lock (_cache/gen_L<L>.lock) makes concurrent calls for the same L wait for the first one (parallel workers
of scripts/paper1_run5.sh).  Resumable: an interrupted run keeps its finished term lists in _cache/terms_L<L>.pkl
and a new call continues from there.
"""
import fcntl
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tb_local  # noqa: E402
from threebody import angular  # noqa: E402
from threebody.coupled import CoupledBasis, block  # noqa: E402


def have(L):
    try:
        d = pickle.load(open(tb_local.LOCAL_CACHE, "rb"))
    except Exception:
        return False
    return any(len(k) > 1 and isinstance(k[1], tuple) and len(k[1]) == L for k in d)


for L in [int(x) for x in sys.argv[1:]]:
    gl = open(os.path.join(os.path.dirname(tb_local.LOCAL_CACHE), f"gen_L{L}.lock"), "w")
    fcntl.flock(gl, fcntl.LOCK_EX)          # released when the process exits
    if have(L):
        print(f"L={L}: already in the local cache", flush=True)
        continue
    t0 = time.perf_counter()
    own = os.path.join(os.path.dirname(tb_local.LOCAL_CACHE), f"terms_L{L}.pkl")
    angular._CACHE_FILE, angular._disk = own, None
    for sym in ("S", "A"):
        CoupledBasis(L, block(2, (1.0, 1.0, 1.0)), sym)
    with open(tb_local.LOCAL_CACHE + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        main = pickle.load(open(tb_local.LOCAL_CACHE, "rb"))
        main.update(pickle.load(open(own, "rb")))
        tmp = tb_local.LOCAL_CACHE + f".{os.getpid()}.tmp"
        pickle.dump(main, open(tmp, "wb"))
        os.replace(tmp, tb_local.LOCAL_CACHE)
    angular._CACHE_FILE, angular._disk = tb_local.LOCAL_CACHE, None
    print(f"L={L}: term lists ready ({time.perf_counter() - t0:.0f} s)", flush=True)
