"""Mass scans over the {X+, X+, e-} family (or any user-defined family)."""
from __future__ import annotations

import csv
import os
import time

import numpy as np

from .convergence import convergence_table, optimize_alphas, default_alphas
from .operators import PAIRS
from .system import symmetric_family


def scan_masses(Ms, Ns, nstates=4, N_opt=12, symmetry="S", warm_start=True, observables=True,
                save_dir=None, system_factory=symmetric_family, alphas=None, grid=None, verbose=True,
                opt_state=0, weights=(1, 1, 1)):
    """For each M: optimize the scale parameters at N_opt (deterministic Nelder-Mead,
    warm-started from the previous M if warm_start), run a convergence table over
    Ns at those fixed scales, extrapolate, and compute observables from the largest basis.

    opt_state selects which root's energy is minimized when choosing the scales
    (0 = ground state; use 1 to describe a weakly bound excited state better).
    Returns a list of row dicts (one per M).
    """
    rows = []
    prev = None
    for M in Ms:
        t0 = time.perf_counter()
        system = system_factory(M)
        if alphas is not None:
            al = tuple(alphas)
            nev = 0
        else:
            x0 = prev if (warm_start and prev is not None) else default_alphas(system)
            al, _, nev = optimize_alphas(system, N_opt, symmetry, state=opt_state, x0=x0, weights=weights)
            prev = al
        tab = convergence_table(system, al, Ns, symmetry, nstates, weights)
        res = tab["last"]
        Eth = system.threshold()
        row = dict(M=M, threshold=Eth, opt_state=opt_state, alpha0=al[0], alpha1=al[1], alpha2=al[2], N_opt=N_opt,
                   N_max=Ns[-1], size=res.basis.size, S_cond=res.info["S_cond_scaled"],
                   n_below_threshold=res.n_below_threshold(),
                   n_below_is_capped_by_nstates=res.n_below_threshold() >= len(res.energies))
        for s in range(nstates):
            E_ext, err, rho, flag = tab["extrap"][s]
            row[f"E{s}"] = res.energies[s] if s < len(res.energies) else np.nan
            row[f"E{s}_ext"] = E_ext
            row[f"E{s}_err"] = err
            row[f"E{s}_flag"] = flag or "ok"
            row[f"E{s}_minus_Eth"] = row[f"E{s}"] - Eth
        row["binding_energy"] = Eth - tab["extrap"][0][0]
        row["dE0_last"] = tab["rows"][-2]["energies"][0] - tab["rows"][-1]["energies"][0] if len(Ns) > 1 else np.nan
        if observables:
            obs = res.observables(0)
            for k, v in obs.items():
                row[f"obs0_{k}"] = v
        row["time_s"] = time.perf_counter() - t0
        row["n_alpha_evals"] = nev
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            res.save(os.path.join(save_dir, f"M_{M:.10g}.npz"))
            if grid is not None:
                write_gr(os.path.join(save_dir, f"gr_M_{M:.10g}.txt"), res, grid)
        rows.append(row)
        if verbose:
            print(f"M={M:<14.10g} E0={row['E0']:.12f} ext={row['E0_ext']:.12f}+/-{row['E0_err']:.1e} "
                  f"Eth={Eth:.10f} n_below={row['n_below_threshold']} ({row['time_s']:.1f}s)", flush=True)
    return rows


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})


def write_gr(path, res, grid, state=0):
    """Pair distribution functions g_ij(r) (int 4 pi r^2 g dr = 1) on a grid, all three pairs."""
    grid = np.asarray(grid, dtype=float)
    cols = [grid] + [res.g_of_r(i, j, grid, state) for (i, j) in PAIRS]
    hdr = "r  g01  g02  g12   (normalized: int 4 pi r^2 g(r) dr = 1); state %d" % state
    np.savetxt(path, np.column_stack(cols), header=hdr)
