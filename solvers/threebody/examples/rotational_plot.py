"""Figure for the rotational tower: reads <outdir>/rotational_table.csv, writes <outdir>/fig_rotational.png."""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..", "..",
                                                         "paper_data", "rotational")
MC1, MC2 = 4.2082827, 11.4894217
rows = list(csv.DictReader(open(os.path.join(OUT, "rotational_table.csv"))))
f = lambda k: np.array([float(r[k]) for r in rows])
yes = lambda k: np.array([r[k] == "True" for r in rows])
M, d1, d2 = f("M"), f("d1"), f("d2")
b1, b2 = yes("bound01"), yes("bound02")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
a = ax[0]
a.semilogx(M[b1], M[b1] * d1[b1], "o-", color="C0", label="exact  M·[E(1)−E(0)]")
a.semilogx(M[b2], M[b2] * d2[b2] / 3, "s-", color="C1", label="exact  M·[E(2)−E(0)]/3")
a.semilogx(M, M * f("BO_d1"), "--", color="C0", alpha=0.6, label="BO  M·ΔE(J=1)")
a.semilogx(M, M * f("BOad_d1"), ":", color="C0", alpha=0.9, label="BO+ad  M·ΔE(J=1)")
a.semilogx(M, M * f("BO_d2") / 3, "--", color="C1", alpha=0.6, label="BO  M·ΔE(J=2)/3")
for m, lab in ((MC1, "L=1 bound"), (MC2, "L=2 bound")):
    a.axvline(m, color="gray", lw=0.8)
    a.text(m * 1.05, 0.6, lab, rotation=90, fontsize=8, color="gray", transform=a.get_xaxis_transform())
a.set_xlabel("M (electron masses)")
a.set_ylabel("M × spacing (hartree × m_e)")
a.set_title("Rotational spacing × M (rigid rotor: 2<1/R²>)")
a.legend(fontsize=7)
a = ax[1]
both = b1 & b2
a.semilogx(M[both], (d2 / d1)[both], "o-", label="exact")
a.semilogx(M, f("BO_ratio"), "--", label="BO")
a.axhline(3, color="k", lw=0.8, label="rigid rotor 3")
a.axvline(MC2, color="gray", lw=0.8)
a.set_xlabel("M")
a.set_ylabel("[E(2)−E(0)] / [E(1)−E(0)]")
a.set_title("Rotational ratio (both levels bound)")
a.legend(fontsize=8)
a = ax[2]
sel = b1 & np.isfinite(f("BOad_d1"))
a.loglog(M[sel], np.abs(d1[sel] / f("BO_d1")[sel] - 1), "o-", label="|d1/d1(BO) − 1|")
a.loglog(M[sel], np.abs(d1[sel] / f("BOad_d1")[sel] - 1), "s-", label="|d1/d1(BO+ad) − 1|")
mm = np.array([M[sel].min(), M[sel].max()])
a.loglog(mm, 1.0 / mm, "k:", label="1/M")
a.set_xlabel("M")
a.set_ylabel("relative deviation of E(L=1)−E(L=0)")
a.set_title("Exact vs BO rotational spacing")
a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_rotational.png"), dpi=130)
print("saved", os.path.join(OUT, "fig_rotational.png"))
