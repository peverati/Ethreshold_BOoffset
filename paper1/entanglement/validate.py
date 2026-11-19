"""Validation of entcore (writes validation.txt).  *Last updated: 2026-11-13*

1. separable Gaussian and separable exponential product: S = 0;
2. coupled Gaussian exp(-a R^2 - b rho^2 - c R.rho): analytic S = 3 S_1D (1D two-mode Gaussian);
   the 1D formula is itself checked by a direct 1D grid SVD;
3. exact state at M = 1: quadrature norm, base vs fine grid.
"""
import math
import sys
import time

import numpy as np

import entcore as ec

out = []


def log(s):
    print(s, flush=True)
    out.append(s)


def gauss_1d(a, b, c):
    A = 2 * (a - c * c / (8 * b))
    B = c * c / (4 * b)
    xi = B / (A + math.sqrt(A * A - B * B))
    S = -math.log(1 - xi) - xi * math.log(xi) / (1 - xi)
    return S, (1 - xi) / (1 + xi)


def gauss_1d_grid(a, b, c, n=400, L=12.0):
    x, w = ec.composite_gl(np.linspace(-L, L, 41), n // 40)
    X, Y = np.meshgrid(x, x, indexing="ij")
    psi = np.exp(-a * X * X - b * Y * Y - c * X * Y)
    A = np.sqrt(w)[:, None] * psi * np.sqrt(w)[None, :]
    lam = np.linalg.svd(A, compute_uv=False) ** 2
    lam /= lam.sum()
    lam = lam[lam > 0]
    return float(-np.sum(lam * np.log(lam))), float(np.sum(lam ** 2))


def make(f):
    def amp(R, r, u):
        rr, uu = np.meshgrid(r, u, indexing="ij")
        return f(R, rr, uu)
    return amp


def run_test(name, f, grid, Sref=None, Pref=None):
    A, direct = ec.partial_waves(make(f), grid)
    sp = ec.spectrum(A, grid.ls)
    rm = ec.radial_measures(A, grid.ls)
    log(f"{name}: S={sp['S']:.12e} S_lin={sp['S_lin']:.12e} S_orient={sp['S_orient']:.6e} S_rad={sp['S_rad']:.6e} "
        f"w0={sp['w'][0]:.12f} E_N(rad)={rm['E_N']:.3e} MI={rm['MI']:.3e} sum|A|^2/direct={sp['norm_sum']/direct.sum():.14f}")
    if Sref is not None:
        log(f"    reference S={Sref:.12e} purity={Pref:.12e}  dS={sp['S']-Sref:.2e}  dP={sp['purity']-Pref:.2e}")
    return sp


if __name__ == "__main__":
    t0 = time.time()
    g = ec.Grid(0.0, 12.0, 12.0, (16, 12), (16, 12), 128, symmetric=True, lmax=40)
    run_test("separable Gaussian exp(-0.6R^2-0.9rho^2)", lambda R, r, u: np.exp(-0.6 * R * R - 0.9 * r * r), g)
    run_test("separable exponential exp(-R-2rho)", lambda R, r, u: np.exp(-R - 2 * r), g)
    a, b, c = 0.5, 0.7, 0.6
    S1, P1 = gauss_1d(a, b, c)
    S1g, P1g = gauss_1d_grid(a, b, c)
    log(f"1D two-mode Gaussian a,b,c={a},{b},{c}: formula S1={S1:.14f} P1={P1:.14f}; grid SVD S1={S1g:.14f} P1={P1g:.14f}")
    gn = ec.Grid(0.0, 10.0, 10.0, (16, 12), (16, 12), 160, symmetric=False, lmax=40)
    # R.rho = R rho u
    run_test("coupled Gaussian exp(-aR^2-b rho^2-c R.rho)",
             lambda R, r, u: np.exp(-a * R * R - b * r * r - c * R * r * u), gn, 3 * S1, P1 ** 3)
    # same in an exchange-symmetric form: cosh(c R.rho) part only has even l -> not analytic; skip
    # exact state, M = 1
    for M in (1.0,):
        res, b_ = ec.load_exact(M)
        amp = ec.exact_amp(res)
        for lev in ("base", "fine"):
            grid = ec.grid_for(M, lev)
            A, direct = ec.partial_waves(amp, grid)
            sp = ec.spectrum(A, grid.ls)
            rm = ec.radial_measures(A, grid.ls)
            log(f"exact M={M} {lev}: direct norm={direct.sum():.12f} sum_l={sp['norm_sum']:.12f} S={sp['S']:.12f} "
                f"S_lin={sp['S_lin']:.12f} S_orient={sp['S_orient']:.10f} S_rad={sp['S_rad']:.10f} "
                f"w0,w2,w4={sp['w'][0]:.10f},{sp['w'][1]:.3e},{sp['w'][2]:.3e} w_lmax={sp['w'][-1]:.2e} "
                f"E_N={rm['E_N']:.8f} MI={rm['MI']:.8f} k={rm['k']},{rm['kp']}")
    log(f"time {time.time()-t0:.1f} s")
    open("validation.txt", "w").write("*Last updated: 2026-11-13*\n\n" + "\n".join(out) + "\n")
