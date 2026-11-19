"""Deterministic choice of scale parameters, convergence tables and extrapolation."""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import minimize

from .basis import PerimetricBasis
from .solver import MatrixSet, solve
from .system import ThreeBodySystem


def default_symmetry(system: ThreeBodySystem):
    return "S" if system.identical is not None else None


def default_alphas(system: ThreeBodySystem):
    """Hydrogenic rule: c_ij = |Z_ij| mu_ij for attractive pairs (0 otherwise),
    alpha_k = c_kj + c_kl.  Exact for separable (one-interaction-off) products."""
    c = {}
    for (i, j), z in system.coupling_dict().items():
        mu = system.reduced_mass(i, j)
        c[(i, j)] = abs(z) * mu if (z < 0 and np.isfinite(mu)) else 0.0
    a = []
    for k in range(3):
        a.append(sum(v for p, v in c.items() if k in p))
    a = np.array(a)
    a[a <= 0] = max(a.max(), 1e-3) * 0.5
    return tuple(a)


def _alphas_from_x(x, symmetric):
    e = np.exp(x)
    return (e[0], e[0], e[1]) if symmetric else tuple(e)


def optimize_alphas(system: ThreeBodySystem, N: int, symmetry="auto", state=0, x0=None,
                    weights=(1, 1, 1), xatol=1e-4, fatol=1e-13, maxiter=400, verbose=False):
    """Minimize the variational energy of `state` w.r.t. the scale parameters at basis size N.

    Deterministic (Nelder-Mead from a deterministic start, in log alpha).
    Returns (alphas, energy, number_of_evaluations).
    """
    if symmetry == "auto":
        symmetry = default_symmetry(system)
    symmetric = symmetry is not None
    a0 = default_alphas(system) if x0 is None else tuple(x0)
    x_init = np.log([a0[0], a0[2]] if symmetric else list(a0))
    count = [0]

    def f(x):
        al = _alphas_from_x(x, symmetric)
        b = PerimetricBasis(N, al, symmetry, weights)
        try:
            e = solve(system, b, nstates=state + 1).energies[state]
        except Exception:  # pragma: no cover - pathological alphas
            return 1e3
        count[0] += 1
        if verbose:
            print(f"  alphas={np.round(al, 6)} E={e:.12f}")
        return e

    res = minimize(f, x_init, method="Nelder-Mead",
                   options=dict(xatol=xatol, fatol=fatol, maxiter=maxiter,
                                initial_simplex=x_init + np.vstack([np.zeros(len(x_init)), 0.2 * np.eye(len(x_init))])))
    return _alphas_from_x(res.x, symmetric), float(res.fun), count[0]


def extrapolate(E):
    """Geometric-tail extrapolation of a monotonically decreasing sequence E[0..n].

    With d_n = E_{n-1} - E_n and rho = d_n / d_{n-1}, the tail is
    d_n rho / (1 - rho).  The error estimate is max(|tail|, |change of the
    extrapolated value from the previous step|): deliberately conservative.
    Returns (E_ext, err, rho, flag) where flag is '' (ok), 'saturated' (last
    change at rounding level), 'nonmonotonic',
    'rho>=1' (not yet in a convergent regime) or 'short'.
    """
    E = np.asarray(E, dtype=float)
    if len(E) < 3:
        return E[-1], (abs(E[-1] - E[-2]) if len(E) == 2 else np.nan), np.nan, "short"
    d = -np.diff(E)  # positive if decreasing
    floor = 100 * np.finfo(float).eps * abs(E[-1])  # rounding noise of the eigenvalues
    if abs(d[-1]) <= floor:
        return E[-1], max(abs(d[-1]), abs(d[-2]) if abs(d[-2]) > floor else floor), 0.0, "saturated"
    if d[-1] < 0 or d[-2] < 0:
        return E[-1], abs(d[-1]) + abs(d[-2]), np.nan, "nonmonotonic"
    if d[-1] == 0:
        return E[-1], 0.0, 0.0, ""
    rho = d[-1] / d[-2]
    if not (0 < rho < 1):
        return E[-1], 10 * d[-1], rho, "rho>=1"
    tail = d[-1] * rho / (1 - rho)
    Eext = E[-1] - tail
    err = tail
    if len(E) >= 4 and d[-3] > 0 and 0 < d[-2] / d[-3] < 1:
        rp = d[-2] / d[-3]
        Eprev = E[-2] - d[-2] * rp / (1 - rp)
        err = max(err, abs(Eext - Eprev))
    return Eext, err, rho, ""


def convergence_table(system: ThreeBodySystem, alphas, Ns, symmetry="auto", nstates=3, weights=(1, 1, 1),
                      keep_last=True):
    """Energies of the lowest nstates for each N (fixed alphas).

    Returns dict(rows=[...], extrap=[(E_ext, err, rho, flag) per state], last=Result).
    """
    if symmetry == "auto":
        symmetry = default_symmetry(system)
    rows = []
    last = None
    for N in Ns:
        t0 = time.perf_counter()
        b = PerimetricBasis(N, alphas, symmetry, weights)
        res = solve(system, b, nstates=nstates)
        rows.append(dict(N=N, size=b.size, energies=res.energies.copy(), time=time.perf_counter() - t0,
                         S_cond=res.info["S_cond_scaled"], dropped=res.info["dropped"]))
        last = res
    ext = []
    for s in range(nstates):
        ext.append(extrapolate([r["energies"][s] for r in rows if len(r["energies"]) > s]))
    return dict(rows=rows, extrap=ext, last=last if keep_last else None, threshold=system.threshold(),
                alphas=tuple(alphas))


def format_table(tab, nstates=None, ref=None):
    rows = tab["rows"]
    ns = len(rows[0]["energies"]) if nstates is None else nstates
    lines = []
    hdr = f"{'N':>3} {'size':>6} " + " ".join(f"{'E'+str(s):>20}" for s in range(ns)) + f" {'dE0':>10} {'t[s]':>7}"
    lines.append(hdr)
    prev = None
    for r in rows:
        e = r["energies"]
        dE = "" if prev is None else f"{prev - e[0]:10.2e}"
        lines.append(f"{r['N']:3d} {r['size']:6d} " + " ".join(f"{x:20.14f}" for x in e[:ns]) + f" {dE:>10} {r['time']:7.2f}")
        prev = e[0]
    lines.append(f"threshold E_th = {tab['threshold']:.14f}; alphas = {np.round(tab['alphas'], 8).tolist()}")
    for s, (E, err, rho, flag) in enumerate(tab["extrap"][:ns]):
        extra = "" if ref is None or s > 0 else f"  (ref {ref:.12f}, diff {E - ref:+.2e})"
        lines.append(f"state {s}: extrapolated {E:.14f} +/- {err:.1e} (rho={rho:.3f}) {flag}{extra}")
    return "\n".join(lines)
