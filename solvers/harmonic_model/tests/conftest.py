import os
import sys

os.environ.setdefault("HARMONIC_BLAS_THREADS", "2")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
