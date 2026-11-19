"""Locations of the smooth features (extrema, inflections in ln M) of the susceptibilities from the
dense windows, with the basis (N vs N-4) and interpolation (spline vs quintic local fit) spread.
Writes features.txt."""
import numpy as np
from scipy.interpolate import CubicSpline
import analyze as an

d = an.load("dense.csv")
win = np.array([k.split("_")[0] for k in d["key"]])
targets = [("w1", "chi_sep", 2, (0.2, 0.3)), ("w1", "chiF", 2, (0.3, 0.38)),
           ("w7", "chi_inv", 2, (1.25, 1.6)), ("w7", "chi_sep", 2, (1.25, 1.6)),
           ("w5", "chi_inv", 1, (5, 9)), ("w5", "chi_sep", 2, (8, 11.5)),
           ("w6", "chi_inv", 1, (15, 40)), ("w6", "chi_sep", 2, (30, 48))]
out = ["window quantity kind  M(N)  M(N-4)  M(local quintic)  value  |deriv order 1,2 at point|"]
for w, k, order, (lo, hi) in targets:
    s = win == w
    L = np.log(d["M"][s])
    res = []
    for col in (k, k + "_N4"):
        cs = CubicSpline(L, d[col][s])
        r = [x for x in cs.derivative(order).roots(extrapolate=False) if np.log(lo) < x < np.log(hi)]
        res.append(np.exp(r[0]) if r else np.nan)
    # local quintic fit around the spline root
    x0 = np.log(res[0])
    sel = np.abs(L - x0) < 0.35
    p = np.polyfit(L[sel] - x0, d[k][s][sel], 5)
    dp = np.polyder(p, order)
    rr = [r.real for r in np.roots(dp) if abs(r.imag) < 1e-12 and abs(r.real) < 0.3]
    mq = np.exp(x0 + min(rr, key=abs)) if rr else np.nan
    cs = CubicSpline(L, d[k][s])
    out.append(f"{w} {k:8s} {'extremum' if order == 1 else 'inflection'}  {res[0]:.4f}  {res[1]:.4f}  {mq:.4f}  "
               f"{float(cs(x0)):.8f}  d1={float(cs(x0, 1)):.3e} d2={float(cs(x0, 2)):.3e}")
open("features.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
