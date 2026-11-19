"""Body-frame (Sanz-Vicario et al. PRA 96, 022503 (2017)) bipartition: scalar |R| versus the
electron (rho, angle to the axis).  Reduced density of |R|: Gamma = sum_l (2l+1) A_l A_l^T / tot,
with A_l the partial-wave radial kernels of entcore.partial_waves (same weights w_l as the lab-frame
R|rho analysis, but the orientational degree of freedom is traced out instead of being a subsystem)."""
import sys, csv, math, os
import numpy as np
import entcore as ec

def body_entropy(M, level):
    res, _ = ec.load_exact(M)
    amp = ec.exact_amp(res)
    grid = ec.grid_for(M, level)
    A, direct = ec.partial_waves(amp, grid)
    mult = 2 * grid.ls + 1
    tot = sum(m * np.sum(a * a) for m, a in zip(mult, A))
    G = sum(m * a @ a.T for m, a in zip(mult, A)) / tot
    lam = np.linalg.eigvalsh(G)[::-1]
    lam = lam[lam > 1e-16]
    S = float(-np.sum(lam * np.log(lam)))
    return dict(M=M, level=level, S_body_nats=S, S_body_bits=S / math.log(2), S_lin=float(1 - np.sum(lam ** 2)),
                lam1=lam[0], lam2=lam[1] if len(lam) > 1 else 0.0, lam3=lam[2] if len(lam) > 2 else 0.0, trace=float(np.sum(lam)))

out = os.path.join(ec.HERE, "results_bodyframe.csv")
done = set()
if os.path.exists(out):
    with open(out) as fh:
        done = {(round(float(r["M"]), 6), r["level"]) for r in csv.DictReader(fh)}
masses = [float(x) for x in sys.argv[1:]] or [1836.15267343]
for M in masses:
    for level in ("base", "fine"):
        if (round(M, 6), level) in done:
            continue
        r = body_entropy(M, level)
        new = not os.path.exists(out)
        with open(out, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(r.keys()))
            if new: w.writeheader()
            w.writerow(r)
        print(r, flush=True)
