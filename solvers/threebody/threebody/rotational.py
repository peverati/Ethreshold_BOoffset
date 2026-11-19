"""Scale optimization, convergence tables and mass scans for L > 0 (natural parity).

Symmetry of the states (exchange P of the identical particles 0, 1; derivation in the
README, section 10): for {X+, X+, e-} the rotational levels of the ground vibronic
manifold (1s sigma_g in the BO picture) have P = (-1)^L, i.e.
L = 0, 2: spatially symmetric ('S'); L = 1: antisymmetric ('A'). Use
`ground_manifold_symmetry(L)`.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import minimize

from .convergence import default_alphas, extrapolate
from .coupled import CoupledBasis, block
from .system import ThreeBodySystem


def ground_manifold_symmetry(L: int) -> str:
    return "S" if L % 2 == 0 else "A"


def optimize_scales(system: ThreeBodySystem, L, N, symmetry, state=0, x0=None, weights=(1, 1, 1),
                    xatol=1e-4, fatol=1e-13, maxiter=300, verbose=False):
    """Nelder-Mead in (log alpha_X, log alpha_e) for one self-mirror total-degree block.

    Deterministic. Returns (alphas, energy, n_evaluations)."""
    a0 = default_alphas(system) if x0 is None else tuple(x0)
    x_init = np.log([a0[0], a0[2]])
    count = [0]

    def f(x):
        al = (np.exp(x[0]), np.exp(x[0]), np.exp(x[1]))
        try:
            b = CoupledBasis(L, block(N, al, weights), symmetry)
            e = b.solve(system, state + 1, vectors=False)[0][state]
        except Exception:  # pragma: no cover
            return 1e3
        count[0] += 1
        if verbose:
            print(f"   alphas={np.round(al, 6)} E={e:.12f}", flush=True)
        return e

    res = minimize(f, x_init, method="Nelder-Mead",
                   options=dict(xatol=xatol, fatol=fatol, maxiter=maxiter,
                                initial_simplex=x_init + np.vstack([np.zeros(2), 0.2 * np.eye(2)])))
    x = res.x
    return (float(np.exp(x[0])), float(np.exp(x[0])), float(np.exp(x[1]))), float(res.fun), count[0]


def convergence_table_L(system, L, symmetry, alphas, Ns, nstates=3, weights=(1, 1, 1), tol=1e-12,
                        observables=False):
    """Energies (fixed scales, nested spaces: monotonic) for increasing N; extrapolation per root."""
    rows = []
    last = None
    for N in Ns:
        t0 = time.perf_counter()
        b = CoupledBasis(L, block(N, alphas, weights), symmetry)
        e, C = b.solve(system, nstates, tol)
        row = dict(N=N, size=b.size, energies=e.copy(), time=time.perf_counter() - t0,
                   dropped=b.project(tol)["dropped"], lam_min=b.project(tol)["lam_min"])
        rows.append(row)
        last = (b, e, C)
    ext = [extrapolate([r["energies"][s] for r in rows]) for s in range(nstates)]
    out = dict(rows=rows, extrap=ext, alphas=tuple(alphas), threshold=system.threshold(), L=L,
               symmetry=symmetry, weights=tuple(weights))
    if observables:
        b, e, C = last
        out["observables"] = observables_L(system, b, C)
    out["last"] = last
    return out


def observables_L(system, b: CoupledBasis, C):
    """<r01>, <r01^2>, <1/r01>, T, V, virial ratio for each column of C (normalized)."""
    im = system.inverse_masses
    z = system.coupling_dict()
    obs = []
    for n in range(C.shape[1]):
        c = C[:, n]
        val = lambda op: float(b.expectation(op, c)[0])
        if b.symmetry is not None:
            T = im[0] * val("T01") + (im[2] * val("T2") if im[2] else 0.0)
            V = z[(0, 1)] * val("V01") + z[(0, 2)] * val("V02+V12")
        else:
            T = sum(im[i] * val(f"T{i}") for i in range(3) if im[i])
            V = sum(zz * val(f"V{i}{j}") for (i, j), zz in z.items() if zz)
        obs.append(dict(r01=val("r01"), r01_2=val("r01^2"), inv_r01=val("1/r01"), T=T, V=V, E=T + V,
                        virial=-V / T))
    return obs


def format_table_L(tab, ref=None):
    lines = [f"L={tab['L']} symmetry={tab['symmetry']} alphas={np.round(tab['alphas'], 8).tolist()} "
             f"weights={list(tab['weights'])} E_th={tab['threshold']:.14f}"]
    prev = None
    for r in tab["rows"]:
        e = r["energies"]
        dE = "" if prev is None else f"{prev - e[0]:10.2e}"
        lines.append(f"{r['N']:3d} {r['size']:6d} " + " ".join(f"{x:20.14f}" for x in e) +
                     f" {dE:>10} {r['time']:7.1f}s drop={r['dropped']}")
        prev = e[0]
    for s, (E, err, rho, flag) in enumerate(tab["extrap"]):
        extra = ""
        if ref is not None and s < len(ref) and ref[s] is not None:
            extra = f"  ref {ref[s]:.15f}  E_last-ref {tab['rows'][-1]['energies'][s] - ref[s]:+.2e}  ext-ref {E - ref[s]:+.2e}"
        lines.append(f"root {s}: ext {E:.14f} +/- {err:.1e} rho={rho:.3f} {flag}{extra}")
    return "\n".join(lines)
