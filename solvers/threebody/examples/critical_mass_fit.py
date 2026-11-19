"""Fits for the critical masses (reads the CSVs written by critical_mass.py).

For each state k:
  1. best binding E_b(M) = E_th - min_basis E_k(M) (a rigorous LOWER bound to the exact
     binding, since each E_k is an upper bound); 'spread' = difference between the best
     and second-best basis of the ladder (convergence indicator).
  2. threshold law, short-range s-wave (ion + neutral polarizable atom, -C4/R^4 tail):
     kappa = sqrt(2 mu E_b) is analytic to leading orders in x = M - M_c, so
         sqrt(E_b) = c1 x + c2 x^2 (+ c3 x^3)
     fitted for (M_c, c_i) on several windows of x; plus a free power law
         E_b = A x^p  (tests p = 2 against e.g. p = 1, which a Coulomb tail would give).
  3. finite-size scaling of the fixed-basis crossing masses M*(alpha_s) -> alpha_s = 0.
  4. size: <r_XX> versus the halo value 1/(2 kappa).
Writes fits_state{k}.txt and fig_state{k}.png into OUTDIR.
"""
import csv
import os
import sys

import numpy as np
from scipy.optimize import curve_fit


def read(path):
    with open(path) as fh:
        return [{k: (float(v) if _isnum(v) else v) for k, v in r.items()} for r in csv.DictReader(fh)]


def _isnum(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


def best_binding(rows):
    Ms = sorted(set(r["M"] for r in rows))
    out = []
    for M in Ms:
        es = sorted(r["E_minus_Eth"] for r in rows if r["M"] == M)
        out.append((M, -es[0], es[1] - es[0]))
    return np.array(out)  # M, Eb, spread


def poly_model(order):
    def f(M, Mc, *c):
        x = M - Mc
        return sum(ci * x ** (i + 1) for i, ci in enumerate(c))
    return f


def fit_sqrt(M, Eb, order, Mc0):
    y = np.sqrt(Eb)
    p0 = [Mc0, y[-1] / (M[-1] - Mc0)] + [0.0] * (order - 1)
    p, cov = curve_fit(poly_model(order), M, y, p0=p0, maxfev=20000)
    res = y - poly_model(order)(M, *p)
    return p, np.sqrt(np.diag(cov)), res


def fit_power(M, Eb, Mc0):
    f = lambda M, Mc, lnA, pw: lnA + pw * np.log(np.clip(M - Mc, 1e-300, None))
    p, cov = curve_fit(f, M, np.log(Eb), p0=[Mc0 - 1e-3, np.log(Eb[-1]) - 2 * np.log(M[-1] - Mc0), 2.0],
                       bounds=([Mc0 - 1.0, -100, 0.5], [M.min() - 1e-12, 100, 4.0]), maxfev=20000)
    return p, np.sqrt(np.diag(cov))


def main(out, states):
    for k in states:
        lines = [f"# Fits for excited symmetric S state k = {k}"]
        er = read(os.path.join(out, f"energies_state{k}.csv"))
        ob = read(os.path.join(out, f"observables_state{k}.csv"))
        cr = read(os.path.join(out, f"crossings_state{k}.csv"))
        data = best_binding(er)
        lines.append("\n## best binding E_b = E_th - min E_k (lower bound to exact binding), spread to 2nd-best basis")
        lines.append(f"{'M':>14} {'E_b':>14} {'spread':>10} {'rel':>8}")
        for M, Eb, sp in data:
            lines.append(f"{M:14.9f} {Eb:14.6e} {sp:10.2e} {sp / abs(Eb):8.1e}")
        # usable: bound and spread < 1% of E_b
        ok = (data[:, 1] > 0) & (data[:, 2] < 1e-2 * data[:, 1])
        M, Eb = data[ok, 0], data[ok, 1]
        std = [r for r in cr if r.get("variant", "standard") == "standard"]
        Mstar = min(r["M_star"] for r in std if np.isfinite(r["M_star"]))
        lines.append(f"\nrigorous upper bound (best fixed-basis crossing): M_c <= {Mstar:.9f}")
        lines.append("\n## threshold-law fits sqrt(E_b) = sum_i c_i (M - M_c)^i on windows (M - M* <= xmax)")
        results = []
        for xmax in (0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
            sel = (M - Mstar) <= xmax + 1e-12
            for order in (1, 2, 3):
                if sel.sum() < order + 2:
                    continue
                try:
                    p, e, res = fit_sqrt(M[sel], Eb[sel], order, Mstar)
                except RuntimeError:
                    continue
                rms = np.sqrt(np.mean(res ** 2)) / np.sqrt(Eb[sel]).max()
                results.append((xmax, order, p[0], e[0], rms))
                lines.append(f"xmax={xmax:5.2f} order={order} n={sel.sum():2d}  M_c={p[0]:.9f} +/- {e[0]:.1e}  "
                             f"c={np.array2string(np.array(p[1:]), precision=6)}  rel.rms={rms:.1e}")
        lines.append("\n## free power law E_b = A (M - M_c)^p")
        for xmax in (0.05, 0.1, 0.3, 1.0):
            sel = (M - Mstar) <= xmax + 1e-12
            if sel.sum() < 4:
                continue
            try:
                p, e = fit_power(M[sel], Eb[sel], Mstar)
                lines.append(f"xmax={xmax:5.2f} n={sel.sum():2d}  M_c={p[0]:.9f} +/- {e[0]:.1e}  p={p[2]:.5f} +/- {e[2]:.1e}")
            except (RuntimeError, ValueError) as exc:
                lines.append(f"xmax={xmax}: fit failed ({exc})")
        lines.append("\n## fixed-basis crossing masses M*(basis) (each a rigorous upper bound)")
        for r in cr:
            lines.append(f"alpha_s={r['alpha_s']:<8g} {r.get('variant', 'standard'):24s} size={int(r['size'])} "
                         f"kept={int(r['kept'])} M*={r['M_star']:.9f}")
        a = np.array([r["alpha_s"] for r in std]); ms = np.array([r["M_star"] for r in std])
        o = np.argsort(a); a, ms = a[o], ms[o]
        if len(a) >= 3:
            for npts in (3, 4):
                pp = np.polyfit(a[:npts], ms[:npts], 1)
                lines.append(f"linear extrapolation of M* to alpha_s -> 0 using {npts} smallest alpha_s: {pp[1]:.9f}")
        lines.append("\n## size of the state (best basis): <r_XX>, rms, halo estimate 1/(2 kappa)")
        lines.append(f"{'M':>14} {'E_b':>12} {'<r_XX>':>12} {'rms r_XX':>12} {'1/(2kappa)':>12} {'ratio':>8} {'<r_Xe>':>10}")
        for r in ob:
            lines.append(f"{r['M']:14.9f} {r['binding']:12.4e} {r['r01']:12.4f} {r['r01_rms']:12.4f} "
                         f"{r['one_over_2kappa']:12.4f} {r['r01'] / r['one_over_2kappa']:8.4f} {r['rXe']:10.4f}")
        txt = "\n".join(lines) + "\n"
        print(txt)
        with open(os.path.join(out, f"fits_state{k}.txt"), "w") as fh:
            fh.write(txt)
        _plot(out, k, data, results, ob, Mstar)


def _plot(out, k, data, results, ob, Mstar):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    pos = data[:, 1] > 0
    ax[0].plot(data[pos, 0], np.sqrt(data[pos, 1]), "o", label="best variational")
    good = [r for r in results if r[1] == 3 and abs(r[0] - 0.1) < 1e-9] or [r for r in results if r[1] == 2]
    if good:
        Mc = good[0][2]
        ax[0].axvline(Mc, color="gray", ls=":", label=f"fit M_c={Mc:.5f}")
    ax[0].axvline(Mstar, color="r", ls="--", lw=0.8, label=f"upper bound {Mstar:.5f}")
    ax[0].set_xlabel("M"); ax[0].set_ylabel("sqrt(E_th - E) [hartree^1/2]"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"state {k}: threshold law")
    Mo = np.array([r["M"] for r in ob]); r01 = np.array([r["r01"] for r in ob])
    h = np.array([r["one_over_2kappa"] for r in ob])
    ax[1].loglog(Mo - Mstar + 1e-9, r01, "o", label="<r_XX>")
    ax[1].loglog(Mo - Mstar + 1e-9, h, "-", label="1/(2 kappa)")
    ax[1].set_xlabel("M - M* (upper bound)"); ax[1].set_ylabel("bohr"); ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"fig_state{k}.png"), dpi=120)


if __name__ == "__main__":
    main(sys.argv[1], [int(s) for s in sys.argv[2:]] or [1, 2, 3])
