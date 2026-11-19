"""Point the threebody angular term cache at this folder's PRIVATE copy (_cache/angular_terms_map.pkl, a copy of
paper1/threshold_map/_cache/angular_terms_map.pkl at commit 8a41b65), so nothing here can modify run 5's cache.
Import this module before any threebody use.  Written 2026-11-17."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers", "threebody"))
from threebody import angular  # noqa: E402

LOCAL_CACHE = os.path.join(HERE, "_cache", "angular_terms_map.pkl")
if not os.path.exists(LOCAL_CACHE):
    raise FileNotFoundError(LOCAL_CACHE)
angular._CACHE_FILE = LOCAL_CACHE
angular._disk = None
