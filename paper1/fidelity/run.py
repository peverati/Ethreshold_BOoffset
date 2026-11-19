"""Resumable driver.
  python run.py grid [M ...]      -> results.csv (mass_scans grid; per-mass optimized fixed basis)
  python run.py dense WINDOW      -> dense.csv  (one fixed-alpha basis per window, fine ln M grid)
Each row is appended and flushed; masses already present (same key) are skipped.
Per row: chi_F (perturbative, exact for the basis), FD check (grid only), chi_inv, chi_sep, the
dilation scalars, the same at N-4 (nested sub-basis), and chi_inv of a pure uniform mass scaling
(exact value 0: a basis-quality indicator for the dilation projection)."""
import csv, os, sys, time
import numpy as np
import fidcore as fc

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.8, 0.9, 1.0, 1.05, 1.1, 1.15,
        1.2, 1.5, 2.0, 3.0, 4.5, 6.0, 9.5256, 12.0, 15.0, 20.0, 29.639, 50.0, 100.0, 206.768283, 500.0,
        1000.0, 1836.15267343, 4000.0, 10000.0]
# dense windows: (name, alpha source mass, M_lo, M_hi, d lnM)
WINDOWS = {"w1": (0.3, 0.17, 0.40, 0.025), "w2": (0.5, 0.33, 0.66, 0.02), "w3": (0.7, 0.5, 0.95, 0.02),
           "w4": (1.0, 0.8, 1.35, 0.025),
           "w5": (6.0, 3.0, 12.5, 0.08), "w6": (20.0, 11.5, 50.0, 0.08),
           "w7": (1.5, 1.2, 2.6, 0.04)}
COLS = ["key", "M", "N", "size", "alpha_X", "alpha_e", "weights", "E0", "dE_dlnM", "chiF", "chiF_fd_h02",
        "chiF_fd_h01", "chiF_fd_rich", "chi_inv", "chi_sep", "DD", "RR", "DR", "pD", "pR", "gram_corr",
        "check_n0", "check_psiD", "check_psiDR", "check_sing_sym", "scal_chi_inv",
        "N4", "E0_N4", "chiF_N4", "chi_inv_N4", "chi_sep_N4", "scal_chi_inv_N4", "time"]


def done_keys(fn):
    if not os.path.exists(fn):
        return set()
    return {r["key"] for r in csv.DictReader(open(fn))}


def append(fn, row):
    new = not os.path.exists(fn)
    with open(fn, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        if new:
            w.writeheader()
        w.writerow(row)
        f.flush(); os.fsync(f.fileno())


def evaluate(P, Nsub, M, fd=False, save=None):
    bs, S, T01, T2, V = P.mats(Nsub)
    E, c, H = fc.ground(S, T01, T2, V, M)
    x, dE = fc.derivative(S, H, E, c, -T01 / M)
    chiF = float(x @ S @ x)
    Cc, Cx = bs.full_tensor(c), bs.full_tensor(x)
    q = fc.dilation_quantities(Cc, Cx, bs)
    r = fc.chis(chiF, q)
    xs, _ = fc.derivative(S, H, E, c, -(T01 / M + T2))
    qs = dict(q)
    qs2 = fc.dilation_quantities(Cc, bs.full_tensor(xs), bs, polar_kw=None)
    rs = fc.chis(float(xs @ S @ xs), qs2)
    out = dict(E0=E, dE_dlnM=dE, chiF=chiF, scal_chi_inv=rs["chi_inv"], q=q, **r)
    if fd:
        vals, rich = fc.chi_fd(S, T01, T2, V, M, c)
        out.update(chiF_fd_h02=vals[0], chiF_fd_h01=vals[1], chiF_fd_rich=rich)
    if save:
        np.savez(save, M=M, N=Nsub, alphas=bs.alphas, weights=bs.weights, c=c, x=x, E=E)
    return out


def run_one(fn, key, M, N, al, wts, P=None, fd=False):
    t0 = time.time()
    if P is None:
        P = fc.Problem(N, al, wts)
    w2 = wts[2]
    a = evaluate(P, N, M, fd=fd, save=os.path.join(fc.DATA, f"dpsi_{key}.npz"))
    b = evaluate(P, N - 4, M)
    q = a["q"]
    row = dict(key=key, M=M, N=N, size=P.sub(N)[0].size, alpha_X=al[0], alpha_e=al[2], weights=str(list(wts)),
               E0=a["E0"], dE_dlnM=a["dE_dlnM"], chiF=a["chiF"], chiF_fd_h02=a.get("chiF_fd_h02", ""),
               chiF_fd_h01=a.get("chiF_fd_h01", ""), chiF_fd_rich=a.get("chiF_fd_rich", ""),
               chi_inv=a["chi_inv"], chi_sep=a["chi_sep"], DD=a["DD"], RR=a["RR"], DR=a["DR"], pD=a["pD"],
               pR=a["pR"], gram_corr=a["gram_corr"], check_n0=q["n0"], check_psiD=q["psi_D"],
               check_psiDR=q["psi_DR"], check_sing_sym=q["RR_sing0"] - q["RR_sing1"],
               scal_chi_inv=a["scal_chi_inv"], N4=N - 4, E0_N4=b["E0"], chiF_N4=b["chiF"],
               chi_inv_N4=b["chi_inv"], chi_sep_N4=b["chi_sep"], scal_chi_inv_N4=b["scal_chi_inv"],
               time=round(time.time() - t0, 1))
    append(fn, row)
    print(f"{key}: M={M:g} chiF={a['chiF']:.8f} inv={a['chi_inv']:.8f} sep={a['chi_sep']:.8f} "
          f"(N-4: {b['chiF']:.8f} {b['chi_inv']:.8f} {b['chi_sep']:.8f}) t={row['time']}s", flush=True)
    return P


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "grid":
        fn = os.path.join(HERE, "results.csv")
        masses = [float(m) for m in sys.argv[2:]] or GRID
        done = done_keys(fn)
        for M in masses:
            key = f"M{M:g}"
            if key in done:
                continue
            N, al, wts, _ = fc.scan_basis(M)
            run_one(fn, key, M, N, al, wts, fd=True)
    elif mode == "dense":
        fn = os.path.join(HERE, "dense.csv")
        for wname in sys.argv[2:]:
            Msrc, lo, hi, step = WINDOWS[wname]
            N, al, wts, _ = fc.scan_basis(Msrc)
            masses = np.exp(np.arange(np.log(lo), np.log(hi) + 1e-9, step))
            done = done_keys(fn)
            P = None
            for M in masses:
                key = f"{wname}_M{M:.6f}"
                if key in done:
                    continue
                P = run_one(fn, key, float(M), N, al, wts, P=P, fd=False)
