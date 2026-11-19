"""Solvable-model chi figure. Reads the harmonic_model CSVs produced by the CLI (see README)."""
import sys, csv, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A, B, out = sys.argv[1], sys.argv[2], "models_chi.csv"
def load(f): return {k: np.array([float(r[k]) for r in csv.DictReader(open(f))]) for k in ("M","chi_fixed","chi_inv","chi_sep")}
a, b = load(A), load(B)
with open(out, "w") as fh:
    fh.write("M,A_chiF,A_chi_inv,A_chi_sep,B_g1_chiF,B_g1_chi_inv,B_g1_chi_sep\n")
    for i in range(len(a["M"])):
        fh.write(",".join("%.12g" % v for v in (a["M"][i], a["chi_fixed"][i], a["chi_inv"][i], a["chi_sep"][i], b["chi_fixed"][i], b["chi_inv"][i], b["chi_sep"][i])) + "\n")
c = {k: np.array([float(r[k]) for r in csv.DictReader(open("../fidelity/results.csv"))]) for k in ("M","chiF","chi_inv","chi_sep")}
o = np.argsort(c["M"]); c = {k: v[o] for k, v in c.items()}
col = {"F": "tab:blue", "inv": "tab:red", "sep": "tab:green"}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
for m, d, ls in (("A", a, "--"), ("B (g=1)", b, "-")):
    for k, key, lab in (("F","chi_fixed","chi_F"),("inv","chi_inv","chi_inv"),("sep","chi_sep","chi_sep")):
        if m == "A" and k == "sep": continue  # identically zero, cannot be drawn on a log axis
        ax[0].plot(d["M"], d[key], ls, color=col[k], label="%s, %s" % (m, lab))
for k, key, lab in (("F","chiF","chi_F"),("inv","chi_inv","chi_inv"),("sep","chi_sep","chi_sep")):
    ax[1].plot(c["M"], c[key], "o-", ms=3, color=col[k], label="Coulomb, " + lab)
ax[1].plot(b["M"], b["chi_sep"], ":", color="gray", lw=1, label="model B (g=1), chi_sep")
for x in ax:
    x.axhline(1/32, color="k", ls=":", lw=1); x.axhline(3/64, color="k", ls="-.", lw=1)
    x.text(0.06, 1/32*1.06, "1/32", fontsize=8); x.text(0.06, 3/64*1.06, "3/64", fontsize=8)
    x.set_xscale("log"); x.set_yscale("log"); x.set_xlim(0.05, 1e4); x.set_ylim(1e-5, 2); x.set_xlabel("M"); x.grid(alpha=.3); x.legend(fontsize=7, loc="lower right")
ax[0].set_ylabel("chi (F = 1 - chi h^2/2, h = d ln M)")
ax[0].set_title("Solvable models A (dashed; chi_sep = 0 identically, not drawn) and B (solid)", fontsize=9)
ax[1].set_title("Exact Coulomb {X+,X+,e-} with model B", fontsize=9)
fig.tight_layout(); fig.savefig("fig_models_chi.png", dpi=150)
