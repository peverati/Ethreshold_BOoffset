"""Radial distribution of the first-order weights at the BO+ad critical points, and WKB (v->inf) limits.
Writes weights_cumulative.csv and wkb_limits.txt."""
import os, sys, numpy as np
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import coupling as C
import pandas as pd
res = pd.read_csv(os.path.join(here, "coupling_results.csv"), comment="#")
edges = [0, 1.5, 3, 5, 8, 15, 40, 1e9]
rows = []
for _, q in res[res.c == 1].iterrows():
    L, lam = int(q.L), q.lam_c
    I, r, u, up = C.weights(lam, L, 1)
    lB = lam * I["B"]
    wr = up**2 / lB; wc = L * (L + 1) * u**2 / r**2 / lB; wB = -C.dU(r) * u**2 / I["B"]; wD = C.dD(r) * u**2 / I["B"]
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        rows.append(dict(L=L, M_c=q.M_c, r_lo=a, r_hi=b, w_rad=np.trapz(wr[m], r[m]), w_cen=np.trapz(wc[m], r[m]),
                         w_B=np.trapz(wB[m], r[m]), adiab_shift_part=np.trapz(wD[m], r[m])))
df = pd.DataFrame(rows); df.to_csv(os.path.join(here, "weights_cumulative.csv"), index=False)
print(df.to_string())
# WKB limit (L=0): dlam_ad = int dD/sqrt(-dU) / int sqrt(-dU) over the allowed region
r = np.geomspace(1e-3, 2e5, 2000001); x = -C.dU(r); m = x > 0
num = np.trapz(C.dD(r[m]) / np.sqrt(x[m]), r[m]); den = np.trapz(np.sqrt(x[m]), r[m])
# radial-weight distribution in WKB: proportional to sqrt(-dU)
cum = [np.trapz(np.sqrt(x[m & (r < b)]), r[m & (r < b)]) / den for b in edges[1:-1]]
txt = f"WKB (v->inf, L=0) first-order adiabatic shift dlam_ad = {num/den:.5f}\nWKB radial-weight cumulative fractions at r < {edges[1:-1]}: {np.round(cum,4).tolist()}\ninner turning point r = {r[m][0]:.4f}\n"
print(txt); open(os.path.join(here, "wkb_limits.txt"), "w").write(txt)
