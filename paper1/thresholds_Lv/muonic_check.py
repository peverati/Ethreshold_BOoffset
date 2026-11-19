"""Direct solves at the physical muonic mass ratios (validation against Frolov 2025, arXiv:2503.17389).

Written 2026-11-16 07:37 . Binding E_th - E (units of the light-particle Hartree,
i.e. muon atomic units for a+ a+ mu-) of
  (J,v) = (1,1): L = 1 'A', root 1, J1v1 bases (scales of J1v1_scales.json, 5 diffuse sets)
  (J,v) = (2,0): L = 2 'S', root 0, same compact scales, 3 diffuse sets
at M = m_d/m_mu and m_t/m_mu as used by Frolov 2025 (17.751673, 26.584936; CODATA 2022 per the notes).
Convergence: compact N = 16, 18, 20. Rows appended to muonic_check.csv (resumable).
Frolov 2025 binding (recomputed from his total energies, arXiv:2503.17389):
  dd-mu (1,1) 3.510171e-4; dd-mu (2,0) 1.537261e-2; tt-mu (2,0) 3.069448e-2.
"""
import csv, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_exact as R
from threebody import symmetric_family, family_threshold
from threebody.coupled import CoupledBasis, block

FROLOV = {("dd", 1, 1): 3.510171e-4, ("dd", 2, 0): 1.537261e-2, ("tt", 2, 0): 3.069448e-2}
MASS = {"dd": 17.751673, "tt": 26.584936}
CASES = [("dd", 1, 1), ("dd", 2, 0), ("tt", 2, 0), ("tt", 1, 1)]
sc = json.load(open(os.path.join(HERE, "J1v1_scales.json")))
mu = sc["M_opt"] / (sc["M_opt"] + 1)
fn = os.path.join(HERE, "muonic_check.csv")
cols = ["ion", "J", "v", "M", "N", "n_diffuse", "size", "E", "E_th", "binding", "frolov", "rel_diff", "seconds"]
done = set()
if os.path.exists(fn):
    done = {(r["ion"], r["J"], r["v"], r["N"]) for r in csv.DictReader(open(fn))}
new = not os.path.exists(fn)
fh = open(fn, "a", newline=""); w = csv.DictWriter(fh, fieldnames=cols)
if new: w.writeheader()
for ion, J, v in CASES:
    M = MASS[ion]
    sym = "A" if J % 2 else "S"
    nd = 5 if J == 1 else 3
    for N in (16, 18, 20):
        if (ion, str(J), str(v), str(N)) in done: continue
        t0 = time.perf_counter()
        blocks = [block(N, al) for al in sc["scales"]]
        for (n_c, n_d, ad) in R.DIFFUSE[:nd]:
            blocks += R.diffuse_pair(mu, ad, n_c, n_d)
        b = CoupledBasis(J, blocks, sym)
        e = b.solve(symmetric_family(M), v + 1, 1e-12, vectors=False)[0]
        Eth = family_threshold(M); Eb = Eth - float(e[v]); fr = FROLOV.get((ion, J, v))
        row = dict(ion=ion, J=J, v=v, M=M, N=N, n_diffuse=nd, size=b.size, E=float(e[v]), E_th=Eth, binding=Eb,
                   frolov=fr if fr else "", rel_diff=(Eb / fr - 1) if fr else "", seconds=round(time.perf_counter() - t0, 1))
        w.writerow(row); fh.flush()
        print(f"{ion}mu ({J},{v}) N={N} size={b.size} E_b={Eb:.9e}" + (f" Frolov {fr:.7e} rel {Eb/fr-1:+.2e}" if fr else "") + f" ({row['seconds']}s)", flush=True)
        del b
