"""Validation of fidcore at one mass (python validate.py M). Prints checks; appends to validation.txt."""
import sys, time
import numpy as np
from scipy.optimize import minimize_scalar
import fidcore as fc

M = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
t0 = time.time()
N, al, wts, Escan = fc.scan_basis(M)
P = fc.Problem(N, al, wts)
bs, S, T01, T2, V = P.mats(N)
E, c, H = fc.ground(S, T01, T2, V, M)
lines = [f"=== M = {M:g}  N = {N} weights {wts} size {len(c)}  E = {E:.14f} (scan {Escan:.14f})  build {time.time()-t0:.1f}s"]
dH = -T01 / M
x, dE = fc.derivative(S, H, E, c, dH)
chiF = x @ S @ x
vals, rich = fc.chi_fd(S, T01, T2, V, M, c)
lines.append(f"chi_F perturbative {chiF:.12f}; FD h=0.02 {vals[0]:.12f} h=0.01 {vals[1]:.12f} Richardson {rich:.12f}; dE/dlnM {dE:.12f}")
Cc, Cx = bs.full_tensor(c), bs.full_tensor(x)
q = fc.dilation_quantities(Cc, Cx, bs)
lines.append(f"quadrature: n0 {q['n0']:.14f} (exact 1) ; nx/n0 {q['nx']/q['n0']:.12f} vs x^T S x {chiF:.12f}; <psi|x> {q['psi_x']:.2e}")
lines.append(f"anti-Hermiticity: <psi|D psi> {q['psi_D']:.2e}  <psi|D_R psi> {q['psi_DR']:.2e}")
lines.append(f"singular pieces (equal by 0<->1 symmetry): {q['RR_sing0']:.12f} {q['RR_sing1']:.12f}")
q2 = dict(q)
for kw in (dict(npan=20, nper=24), dict(extra=20)):
    a = fc.polar_singular(Cc, bs, 0, **kw)
    lines.append(f"polar convergence {kw}: {a:.12f} (diff {a-q['RR_sing0']:.2e})")
qq = fc.dilation_quantities(Cc, Cx, bs, extra=16)
lines.append(f"product-grid extra=16: DD {qq['DD']:.12f} vs {q['DD']:.12f}; pR {qq['pR']:.12f} vs {q['pR']:.12f}")
r = fc.chis(chiF, q)
lines.append("chis: " + " ".join(f"{k} {v:.10f}" for k, v in r.items()))
# D norm against exact dilated overlaps
fdD = []
for sg in (0.01, 0.005):
    o = 0.5 * (fc.dilated_overlap(Cc, al, Cc, al, sg) + fc.dilated_overlap(Cc, al, Cc, al, -sg))
    fdD.append(2 * (1 - o) / sg ** 2)
lines.append(f"||D psi||^2: quadrature {r['DD']:.10f}; from exact dilated overlaps {fdD[0]:.10f} {fdD[1]:.10f} Rich {(4*fdD[1]-fdD[0])/3:.10f}")
# uniform mass scaling (pure dilation): all masses * e^h
xs, _ = fc.derivative(S, H, E, c, -(T01 / M + T2))
qs = fc.dilation_quantities(Cc, bs.full_tensor(xs), bs)
chis_s = xs @ S @ xs
rs = fc.chis(chis_s, qs)
lines.append(f"uniform-mass scaling: chi {chis_s:.10f} vs ||Dpsi||^2 {rs['DD']:.10f}; <x|Dpsi>/DD {rs['pD']/rs['DD']:.10f} (exact 1); chi_inv {rs['chi_inv']:.2e} chi_sep {rs['chi_sep']:.2e} (exact 0)")
# finite-separation dilation-invariant fidelity
for h in (0.02, 0.01):
    Fs = []
    for sgn in (1, -1):
        _, c2, _ = fc.ground(S, T01, T2, V, M * np.exp(sgn * h))
        C2 = bs.full_tensor(c2)
        res = minimize_scalar(lambda sg: -abs(fc.dilated_overlap(Cc, al, C2, al, sg)), bounds=(-0.05, 0.05), method="bounded", options=dict(xatol=1e-10))
        Fs.append(-res.fun)
    lines.append(f"finite-separation chi_inv h={h}: {2*(1-np.mean(Fs))/h**2:.10f} (vs {r['chi_inv']:.10f})")
lines.append(f"time {time.time()-t0:.1f}s")
print("\n".join(lines))
open("validation.txt", "a").write("\n".join(lines) + "\n")
