"""Crossover width of the atom-to-molecule change of the {X+,X+,e-} ground state (analysis of existing data).

Reads ../mass_scans/descriptors.csv (Lindemann ratio of the exact ground state, 37 masses 0.05..1e4).
Two definitions of the width in ln M:
 (a) local exponent f(M) = -4 dlnL/dlnM, which goes from 0 (Lindemann maximum, atom-like side)
     to 1 (BO limit L ~ M^-1/4); width = ln M between f = 0.1 and f = 0.9.
 (b) Lindemann ratio from 90% to 10% of its maximum (10% point from the M^-1/4 tail of the last point).
Then the spread of the ground-state markers in units of each width.  Writes results.txt.
"""
import numpy as np, pandas as pd, os
H = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(H, "..", "mass_scans", "descriptors.csv")).sort_values("M").drop_duplicates("M")
M, L = d.M.values, d.lindemann.values
x = np.log(M)

def cross(level, xs, ys):
    i = np.where((ys[:-1] - level) * (ys[1:] - level) <= 0)[0]
    return [float(np.exp(xs[j] + (level - ys[j]) * (xs[j + 1] - xs[j]) / (ys[j + 1] - ys[j]))) for j in i]

xm = (x[1:] + x[:-1]) / 2
f = -4 * np.diff(np.log(L)) / np.diff(x)
fa = {lev: cross(lev, xm, f)[0] for lev in (0.1, 0.5, 0.9)}
wa = np.log(fa[0.9] / fa[0.1])
Lmax = L.max()
b90 = cross(0.9 * Lmax, x, L)[0]
b50 = cross(0.5 * Lmax, x, L)[0]
b10 = (L[-1] * M[-1] ** 0.25 / (0.1 * Lmax)) ** 4
wb = np.log(b10 / b90)
markers = {"ground-state markers, all (Lindemann max .. X kinetic max)": (0.2, 1.11),
           "literature-proposed only (Matyus 0.4-0.8, King 0.4924, Kais-Shi 0.543, Salas 0.98)": (0.4, 0.98),
           "with the sharp thresholds (L=1 4.21 .. v=4 98.27)": (0.2, 98.2655)}
out = ["# Crossover width (ground-state Lindemann ratio), from ../mass_scans/descriptors.csv", "",
       f"(a) local exponent f=-4 dlnL/dlnM: f=0.1 at M={fa[0.1]:.3f}, f=0.5 at M={fa[0.5]:.3f}, f=0.9 at M={fa[0.9]:.2f}; "
       f"width {wa:.2f} in ln M ({wa/np.log(10):.2f} decades)",
       f"(b) Lindemann 90% of max ({Lmax:.4f} at M={M[L.argmax()]:g}) at M={b90:.3f}, 50% at M={b50:.2f}, "
       f"10% at M={b10:.3g} (M^-1/4 tail extrapolation); width {wb:.2f} in ln M ({wb/np.log(10):.2f} decades)", ""]
for k, (lo, hi) in markers.items():
    s = np.log(hi / lo)
    out.append(f"{k}: spread {s:.2f} in ln M = {s/wa:.2f} x width(a) = {s/wb:.2f} x width(b)")
out += ["", "Grid: 37 masses, finite differences between neighbours (f is evaluated at geometric midpoints);",
        "interpolation linear in ln M. The f levels move by a few percent with the grid; the conclusion does not."]
open(os.path.join(H, "results.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
