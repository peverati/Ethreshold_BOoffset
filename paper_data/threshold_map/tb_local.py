"""Local import shim for solvers/threebody (threshold_map; written 2026-11-16).

solvers/ must not be modified by this project folder, but threebody.angular persists its sympy-derived term
lists in solvers/threebody/threebody/_cache/angular_terms_v1.pkl, which holds L <= 4 only.  The L = 5, 6, 7
term lists needed here are written to a LOCAL cache instead: paper_data/threshold_map/_cache/angular_terms_map.pkl,
seeded (once) with a copy of the solver cache.  Import this module before any threebody computation:
    import tb_local  # noqa: F401
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
from threebody import angular  # noqa: E402

LOCAL_CACHE = os.path.join(HERE, "_cache", "angular_terms_map.pkl")
if not os.path.exists(LOCAL_CACHE):
    os.makedirs(os.path.dirname(LOCAL_CACHE), exist_ok=True)
    if os.path.exists(angular._CACHE_FILE):
        shutil.copy2(angular._CACHE_FILE, LOCAL_CACHE)
angular._CACHE_FILE = LOCAL_CACHE      # the disk-cache wrapper reads this module global at call time
angular._disk = None
