"""Non-BO critical masses M_c1..3 of {X+, X+, e-} (symmetric S states) with the multi-scale basis.

Stage "compute": for each excited root k = 1, 2, 3
  * basis family  B(alpha_s) = compact blocks (fixed scales) + D(4, 40, alpha_s) + D(4, 30, 0.03),
    alpha_s on a ladder (more and more diffuse);
  * E_k(M) - E_th(M) on a mass grid for every basis (best = minimum: every value is a
    rigorous upper bound);
  * crossing mass M*(basis) where E_k = E_th at fixed basis: rigorous UPPER bound to M_ck;
  * convergence checks: larger basis (compact N 20, D(6,60) + D(6,40)), orthogonalization tol;
  * <r_XX>, <r_XX^2>, <r_Xe> of the state along M (best basis).
Writes CSV files to OUTDIR.  Stage "fit" is in critical_mass_fit.py.

Run from solvers/threebody:  THREEBODY_BLAS_THREADS=3 python examples/critical_mass.py OUTDIR [states]
"""
import csv
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np  # noqa: E402

from threebody import family_threshold  # noqa: E402
from threebody.threshold import MultiBlockBasis, compact_block, diffuse_block  # noqa: E402

OPS = ("S", "TX", "Te", "V01", "VXe", "r01", "r01^2", "rXe")

# compact scales (alpha_X, alpha_e): ground-state optimum, excited-root optimum (package
# optimize_alphas at N = 12), and for k >= 2 an intermediate set.  M_ref fixes the atom scale
# of the diffuse blocks (mu_a = M_ref/(M_ref+1)).  The basis does not depend on M.
CONFIG = {
    1: dict(M_ref=10.0, compact=[(0.454, 0.98), (1.458, 1.588)], M_lo=9.0, M_hi=10.0),
    2: dict(M_ref=30.0, compact=[(2.554, 1.442), (0.97, 1.2), (0.369, 0.985)], M_lo=28.0, M_hi=30.0),
    3: dict(M_ref=60.0, compact=[(3.487, 1.420), (1.38, 1.2), (0.545, 1.003)], M_lo=56.0, M_hi=60.0),
}
LADDER = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5]
DELTAS = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]


def basis(k, alpha_s, Nc=16, nc=4, nd=40, nd2=30):
    c = CONFIG[k]
    mu = c["M_ref"] / (c["M_ref"] + 1)
    blocks = [compact_block(Nc, ax, ae) for ax, ae in c["compact"]]
    blocks += [diffuse_block(mu, alpha_s, nc, nd), diffuse_block(mu, 0.03, nc, nd2)]
    return MultiBlockBasis(blocks, ops=OPS)


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({a: (repr(v) if isinstance(v, float) else v) for a, v in r.items()})


def run_state(k, out, log):
    c = CONFIG[k]
    t_state = time.perf_counter()
    bases = {}
    cross = []
    for a in LADDER:
        t0 = time.perf_counter()
        mb = basis(k, a)
        pr = mb.project(1e-12)
        try:
            Ms = mb.critical_mass(c["M_lo"], c["M_hi"], state=k)
        except ValueError as exc:
            Ms = float("nan")
            log(f"  k={k} alpha_s={a}: {exc}")
        bases[a] = mb
        cross.append(dict(state=k, alpha_s=a, size=mb.size, kept=pr["kept"], dropped=pr["dropped"],
                          M_star=Ms, build_s=mb.build_time, total_s=time.perf_counter() - t0))
        log(f"  k={k} alpha_s={a:g} size={mb.size} kept={pr['kept']} M*={Ms:.9f} ({time.perf_counter() - t0:.0f}s)")
    Mstar = np.nanmin([r["M_star"] for r in cross])
    grid = [Mstar + d for d in DELTAS]
    rows = []
    obs = []
    for M in grid:
        Eth = family_threshold(M)
        best = None
        for a, mb in bases.items():
            e = mb.energies(M, k + 1)[k]
            rows.append(dict(state=k, M=M, alpha_s=a, E=e, E_minus_Eth=e - Eth))
            if best is None or e < best[1]:
                best = (a, e)
        ob = bases[best[0]].observables(M, state=k)
        mu_ion = M * (M + 1) / (2 * M + 1)
        Eb = Eth - ob["E"]
        kappa = np.sqrt(2 * mu_ion * Eb) if Eb > 0 else float("nan")
        obs.append(dict(state=k, M=M, best_alpha_s=best[0], E=ob["E"], E_th=Eth, binding=Eb, kappa=kappa,
                        r01=ob["r01"], r01_rms=np.sqrt(ob["r01_2"]), rXe=ob["rXe"], one_over_2kappa=1 / (2 * kappa),
                        dE_dM=ob["dE_dM"], dEth_dM=-1 / (2 * (M + 1) ** 2)))
    log(f"  k={k} mass grid done ({time.perf_counter() - t_state:.0f}s so far)")

    # ---- convergence checks at three masses
    conv = []
    checkM = [Mstar + 0.01, Mstar + 0.1, Mstar + 0.5]
    best_alpha = {r["M"]: r["best_alpha_s"] for r in obs}
    for M in checkM:
        a = min(LADDER, key=lambda x: min(bases[x].energies(M, k + 1)[k], 1))
        e = bases[a].energies(M, k + 1)[k]
        for tol in (1e-10, 1e-14):
            et = bases[a].energies(M, k + 1, tol=tol)[k]
            conv.append(dict(state=k, M=M, variant=f"alpha_s={a:g} tol={tol:g}", size=bases[a].size,
                             E_minus_Eth=et - family_threshold(M), diff_from_ref=et - e))
        conv.append(dict(state=k, M=M, variant=f"alpha_s={a:g} reference (tol 1e-12)", size=bases[a].size,
                         E_minus_Eth=e - family_threshold(M), diff_from_ref=0.0))
    for label, kw in [("compact N=20", dict(Nc=20)), ("diffuse nc=6", dict(nc=6)),
                      ("diffuse nd=60, nd2=45", dict(nd=60, nd2=45))]:
        t0 = time.perf_counter()
        for a in (1e-3, 1e-4):
            mb = basis(k, a, **kw)
            for M in checkM:
                e = mb.energies(M, k + 1)[k]
                ref = bases[a].energies(M, k + 1)[k]
                conv.append(dict(state=k, M=M, variant=f"{label}, alpha_s={a:g}", size=mb.size,
                                 E_minus_Eth=e - family_threshold(M), diff_from_ref=e - ref))
            if label == "compact N=20" or a == 1e-4:
                try:
                    Ms = mb.critical_mass(c["M_lo"], c["M_hi"], state=k)
                except ValueError:
                    Ms = float("nan")
                cross.append(dict(state=k, alpha_s=a, size=mb.size, kept=mb.project()["kept"],
                                  dropped=mb.project()["dropped"], M_star=Ms, build_s=mb.build_time,
                                  total_s=time.perf_counter() - t0, variant=label))
        log(f"  k={k} check '{label}' done ({time.perf_counter() - t0:.0f}s)")
    for r in cross:
        r.setdefault("variant", "standard")
    del best_alpha
    write_csv(os.path.join(out, f"energies_state{k}.csv"), rows)
    write_csv(os.path.join(out, f"observables_state{k}.csv"), obs)
    write_csv(os.path.join(out, f"crossings_state{k}.csv"), cross)
    write_csv(os.path.join(out, f"convergence_state{k}.csv"), conv)
    with open(os.path.join(out, f"basis_state{k}.json"), "w") as fh:
        json.dump(dict(config=CONFIG[k], ladder=LADDER, blocks=basis_desc(bases[LADDER[0]]),
                       wall_time_s=time.perf_counter() - t_state), fh, indent=1)
    log(f"state {k} finished in {time.perf_counter() - t_state:.0f}s")


def basis_desc(mb):
    return mb.describe()


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    states = [int(s) for s in sys.argv[2:]] or [1, 2, 3]
    os.makedirs(out, exist_ok=True)
    logf = open(os.path.join(out, "run_log.txt"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"=== critical_mass.py states={states} {time.strftime('%Y-%m-%d %H:%M:%S')}")
    for k in states:
        run_state(k, out, log)


if __name__ == "__main__":
    main()
