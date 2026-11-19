"""Fixed-exponent fits of y = M^1.5 (E_exact - E_BOad) = c6 + c8 k^2 + c10 k^4 (k = M^-1/4), and of
m_v_eff = -(E_exact - E_BOad) M / T_rad^BOad = m0 + m1 k^2 + m2 k^4. Input: kappa6_data.csv."""
import pandas as pd, numpy as np, itertools, os
here = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(os.path.join(here, "kappa6_data.csv"))
out = []
def fit(col, nterm, Mmin, Mmax=1e9):
    s = d[(d.M >= Mmin) & (d.M <= Mmax)]
    k2 = s.kappa.values**2
    A = np.vstack([k2**j for j in range(nterm)]).T
    w = 1 / (s.E_ex_err.values * s.M.values**1.5 + 1e-9)   # weights from exact-energy error
    c, *_ = np.linalg.lstsq(A * w[:, None], s[col].values * w, rcond=None)
    res = s[col].values - A @ c
    return c, len(s), np.max(np.abs(res))
for col in ("y6", "mv_eff"):
    for nterm in (2, 3):
        for Mmin in (100, 206.768283, 500, 1000, 1836.152673):
            for Mmax in (1e4, 4000):
                c, n, r = fit(col, nterm, Mmin, Mmax)
                if n < nterm + 1: continue
                out.append((col, nterm, Mmin, Mmax, n, *c[:2], c[2] if nterm > 2 else np.nan, r))
o = pd.DataFrame(out, columns=["quantity", "nterms", "Mmin", "Mmax", "npts", "c0", "c1", "c2", "maxres"])
o.to_csv(os.path.join(here, "kappa6_fits.csv"), index=False)
print(o.to_string())
for col in ("y6", "mv_eff"):
    s = o[o.quantity == col]
    print(col, "c0 range", s.c0.min(), s.c0.max(), " 3-term:", s[s.nterms == 3].c0.min(), s[s.nterms == 3].c0.max())
