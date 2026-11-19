"""Unattended convergence of the M = 1e4 levels (L = 0, 1, 2) to larger basis size N.

Last updated: 2026-11-14 17:40 (see README.md).

Continues the nested fixed-scale sequences of paper1/mass_scans (same alphas, read from
exact_L0.csv / exact_rot.csv) to larger N. ONE ROW PER (job, N) is appended to results.csv
right after each diagonalization, so a killed run loses at most the N in progress; rerunning
the same command skips every (job, N) already present.  Jobs run cheapest first.

usage:  python run.py            # all jobs
        python run.py --dry      # print the plan and the estimated matrix sizes, run nothing

Run 3 (2026-11-14): L2 N = 46 (size 14400) added.  The standard path (CoupledBasis keeps 5 raw
operator matrices + projected copies + LAPACK work, measured peak ~13 n^2 doubles, i.e. ~22 GB at
n = 14400) would swap on a 16 GB machine, so (job, N) pairs in LEAN use solve_lean(): S is
diagonalized alone, H is assembled directly as ONE matrix (no per-operator storage), projected
and diagonalized; peak ~4 n^2 doubles (~7 GB at n = 14400).  Same function space, same tol,
same energy (checked against the standard path at N = 26 and 30: differences 6e-14 and 2e-13,
rounding level). Measured peak footprint at N = 30: 0.53 GB lean vs 2.03 GB standard.
"""
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SCANS = os.path.join(HERE, "..", "mass_scans")
sys.path.insert(0, SCANS)
os.environ.setdefault("THREEBODY_BLAS_THREADS", "8")
import common  # noqa: E402,F401  (sets sys.path for threebody)
from threebody import symmetric_family, PerimetricBasis, MatrixSet, solve, extrapolate  # noqa: E402
from threebody.rotational import convergence_table_L  # noqa: E402

M = 10000.0
# (name, L, sym, new Ns). Old Ns/energies come from the mass_scans CSVs.
JOBS = [("L0", 0, "S", [52, 56, 60]),
        ("L1", 1, "A", [42, 46]),
        ("L2", 2, "S", [38, 42]),
        # tier 2 (largest matrices, ~12-13k; may swap on a 16 GB machine; runs last so tier 1 is safe)
        ("L0", 0, "S", [64]),
        ("L1", 1, "A", [50]),
        # run 3: one more L = 2 size (14400 functions), memory-lean path
        ("L2", 2, "S", [46])]
LEAN = {("L2", 46)}
OUT = os.path.join(HERE, "results.csv")
COLS = ["job", "M", "L", "sym", "N", "size", "E", "alpha_X", "alpha_e", "weights", "time_s", "finished"]


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
    print(line, flush=True)
    with open(os.path.join(HERE, "run_log.txt"), "a") as fh:
        fh.write(line + "\n")


def previous(L):
    """alphas and old (N, E) sequence of the mass_scans row at M = 1e4."""
    if L == 0:
        rows = [r for r in csv.DictReader(open(os.path.join(SCANS, "exact_L0.csv"))) if float(r["M"]) == M]
        r = rows[-1]
        seq = json.loads(r["E0_seq"])
        Ns = [36, 40, 44, 48]
    else:
        rows = [r for r in csv.DictReader(open(os.path.join(SCANS, "exact_rot.csv")))
                if float(r["M"]) == M and int(r["L"]) == L]
        r = rows[-1]
        seq = json.loads(r["E_seq"])
        Ns = list(range(int(r["N_last"]) - 4 * (len(seq) - 1), int(r["N_last"]) + 1, 4))
    al = (float(r["alpha_X"]), float(r["alpha_X"]), float(r["alpha_e"]))
    return al, tuple(json.loads(r["weights"])), list(zip(Ns, seq))


def solve_lean(system, L, sym, alphas, N, weights, tol=1e-12):
    """Lowest eigenvalue for CoupledBasis(L, block(N), sym) with ~4 n^2 peak memory (see docstring)."""
    import gc
    import numpy as np
    from scipy.linalg import eigh
    from threebody import angular
    from threebody.coupled import CoupledBasis, block
    from threebody.operators import combine
    from threebody.threshold import raw_matrix
    b = CoupledBasis(L, block(N, alphas, weights), sym, ops=("S",))   # assembles S only
    n = b.size
    S = b.mats.pop("S")
    dd = 1.0 / np.sqrt(np.diag(S).copy())
    S *= dd[:, None]
    S *= dd[None, :]
    lam, U = eigh(S, overwrite_a=True, check_finite=False)
    del S
    gc.collect()
    i0 = int(np.searchsorted(lam, tol * lam[-1], side="right"))     # lam ascending: keep lam[i0:]
    X = U[:, i0:]
    X *= dd[:, None]
    X /= np.sqrt(lam[i0:])[None, :]
    dropped = i0
    im = system.inverse_masses
    z = system.coupling_dict()
    assert system.identical is not None and sym is not None

    def hterms(c, d):
        parts = [(im[0], angular.operator_terms("T01", c, d)),
                 (z[(0, 1)], angular.operator_terms("V01", c, d)),
                 (z[(0, 2)], angular.operator_terms("V02+V12", c, d))]
        if im[2] != 0.0:
            parts.append((im[2], angular.operator_terms("T2", c, d)))
        return combine(*parts)

    H = np.zeros((n, n))
    cache = {}
    off = b.offsets
    for I, (c, B, nb) in enumerate(b.units):
        for J in range(I, len(b.units)):
            d, C, nc = b.units[J]
            blk = raw_matrix(hterms(c, d), B, C, cache)
            blk = blk + b.eps * raw_matrix(hterms(c, angular.swap_component(d)), B, C.mirror(), cache)
            blk *= 2.0 * np.outer(nb, nc)
            if J == I:
                blk = 0.5 * (blk + blk.T)
                H[off[I]:off[I + 1], off[I]:off[I + 1]] = blk
            else:
                H[off[I]:off[I + 1], off[J]:off[J + 1]] = blk
                H[off[J]:off[J + 1], off[I]:off[I + 1]] = blk.T
            del blk
    del cache
    T = H @ X
    del H
    gc.collect()
    Hp = X.T @ T
    del T, X, U
    gc.collect()
    Hp = 0.5 * (Hp + Hp.T)
    e = eigh(Hp, eigvals_only=True, subset_by_index=[0, 0], overwrite_a=True, check_finite=False)
    return float(e[0]), n, dropped


def have():
    if not os.path.exists(OUT):
        return {}
    out = {}
    for r in csv.DictReader(open(OUT)):
        out[(r["job"], int(r["N"]))] = float(r["E"])
    return out


def main(dry=False):
    s = symmetric_family(M)
    done = have()
    new = not os.path.exists(OUT)
    fh = open(OUT, "a", newline="")
    wr = csv.DictWriter(fh, fieldnames=COLS)
    if new:
        wr.writeheader()
    for name, L, sym, Ns in JOBS:
        al, w, old = previous(L)
        log(f"{name}: alphas={al} weights={w} old N,E = {old}")
        for N in Ns:
            if (name, N) in done:
                log(f"{name} N={N}: already done, E={done[(name, N)]:.14f}")
                continue
            if dry:
                log(f"{name} N={N}: TO DO" + (" (memory-lean path)" if (name, N) in LEAN else ""))
                continue
            t0 = time.time()
            if (name, N) in LEAN:
                E, size, _ = solve_lean(s, L, sym, al, N, w)
            elif L == 0:
                b = PerimetricBasis(N, al, "S", w)
                r = solve(s, b, nstates=1, matrices=MatrixSet(b))
                E, size = float(r.energies[0]), b.size
                del r, b
            else:
                tab = convergence_table_L(s, L, sym, al, [N], nstates=1, weights=w, observables=False)
                E, size = float(tab["rows"][0]["energies"][0]), tab["rows"][0]["size"]
            wr.writerow(dict(job=name, M=M, L=L, sym=sym, N=N, size=size, E=f"{E:.16g}", alpha_X=al[0],
                             alpha_e=al[2], weights=json.dumps(w), time_s=round(time.time() - t0, 1),
                             finished=time.strftime("%Y-%m-%d %H:%M")))
            fh.flush()
            done[(name, N)] = E
            log(f"{name} N={N} size={size}: E={E:.14f} t={time.time() - t0:.0f}s")
        seq = [e for _, e in old] + [done[k] for k in sorted(k for k in done if k[0] == name)]
        if len(seq) >= 3:
            e = extrapolate(seq)
            log(f"{name} SUMMARY: E_last={seq[-1]:.14f} ext={e[0]:.14f} +/- {e[1]:.1e} rho={e[2]:.3f} {e[3]}")
    log("ALL JOBS FINISHED" if not dry else "dry run done")


if __name__ == "__main__":
    main(dry="--dry" in sys.argv)
