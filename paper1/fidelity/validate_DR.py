"""Independent brute-force check of the separate-dilation scalars: direct quadrature in Jacobi
coordinates (R, rho, theta), D_R psi = d/de psi(R e^e, rho, theta)|_0 + 3/2 psi by central FD.
Usage: python validate_DR.py M  (compares with fidcore polar/product results)."""
import sys
import numpy as np
from scipy.special import roots_legendre
import fidcore as fc
from threebody import laguerre

M = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
N, al, wts, _ = fc.scan_basis(M)
P = fc.Problem(N, al, wts)
bs, S, T01, T2, V = P.mats(N)
E, c, H = fc.ground(S, T01, T2, V, M)
x, _ = fc.derivative(S, H, E, c, -T01 / M)
Cc, Cx = bs.full_tensor(c), bs.full_tensor(x)
q = fc.dilation_quantities(Cc, Cx, bs)
ref = dict(RR=q["RR"] / q["n0"], pR=q["pR"] / q["n0"], DR=q["DR"] / q["n0"])

def psi(C, r01, r02, r12):
    s0 = r01 + r02 - r12; s1 = r01 + r12 - r02; s2 = r02 + r12 - r01
    F0 = laguerre.evaluate(bs.nmax[0], al[0], s0); F1 = laguerre.evaluate(bs.nmax[1], al[1], s1)
    F2 = laguerre.evaluate(bs.nmax[2], al[2], s2)
    T = np.einsum("abc,cq->abq", C, F2, optimize=True)
    T = np.einsum("abq,bq->aq", T, F1, optimize=True)
    return np.einsum("aq,aq->q", T, F0, optimize=True)

def comp_gl(edges, n):
    x, w = roots_legendre(n); X, W = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        X.append(0.5*(hi-lo)*x+0.5*(hi+lo)); W.append(0.5*(hi-lo)*w)
    return np.concatenate(X), np.concatenate(W)

L = 1.0 / min(al)   # length scale
Rn, Rw = comp_gl(np.concatenate([[0], np.geomspace(1e-3*L, 60*L, 24)]), 12)
th, thw = comp_gl(np.linspace(0, np.pi, 5), 16)
ct, st = np.cos(th), np.sin(th)
acc = dict(n0=0.0, RR=0.0, pR=0.0, DR=0.0)
eps = 1e-4
for R, wR in zip(Rn, Rw):
    edges = np.unique(np.concatenate([[0], np.geomspace(1e-3*L, 60*L, 20), [R/2]]))
    edges = np.sort(edges)
    rho, rw = comp_gl(edges, 10)
    RH, CT = np.meshgrid(rho, ct, indexing="ij"); W = np.outer(rw * rho**2, thw * st) * wR * R**2
    def dist(Rv):
        r02 = np.sqrt(np.maximum(RH**2 + Rv**2/4 + RH*Rv*CT, 0)); r12 = np.sqrt(np.maximum(RH**2 + Rv**2/4 - RH*Rv*CT, 0))
        return np.full_like(RH, Rv).ravel(), r02.ravel(), r12.ravel()
    p0 = psi(Cc, *dist(R)); px = psi(Cx, *dist(R))
    pp = psi(Cc, *dist(R*np.exp(eps))); pm = psi(Cc, *dist(R*np.exp(-eps)))
    p2p = psi(Cc, *dist(R*np.exp(2*eps))); p2m = psi(Cc, *dist(R*np.exp(-2*eps)))
    DRp = (8*(pp-pm) - (p2p-p2m)) / (12*eps) + 1.5*p0
    # D psi: scale R and rho together = Euler in distances: use psi(lambda r)
    Dp = (8*(psi(Cc, *[d*np.exp(eps) for d in dist(R)]) - psi(Cc, *[d*np.exp(-eps) for d in dist(R)]))
          - (psi(Cc, *[d*np.exp(2*eps) for d in dist(R)]) - psi(Cc, *[d*np.exp(-2*eps) for d in dist(R)]))) / (12*eps) + 3*p0
    Wf = W.ravel()
    acc["n0"] += Wf @ p0**2; acc["RR"] += Wf @ DRp**2; acc["pR"] += Wf @ (px*DRp); acc["DR"] += Wf @ (Dp*DRp)
out = {k: acc[k]/acc["n0"] for k in ("RR", "pR", "DR")}
line = f"M={M:g} brute-force Jacobi: " + " ".join(f"{k} {out[k]:.9f} (fidcore {ref[k]:.9f}, rel {out[k]/ref[k]-1:.1e})" for k in out)
print(line); open("validation.txt", "a").write(line + "\n")
