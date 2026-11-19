"""Zero-energy (threshold) coupling eigenproblems for {A+ (M), B+ (rM), e-}: BO, BO+ad, two-channel g/u model.

Written 2026-11-16.  Theory: THEORY.md.  Scaled units r = mu_e R, lam = M_eff/mu_e (linear coupling constant),
M = (1+r) lam/(2r) - 1/(1+r) (bo_h2plus.unequal.M_of_lam).  Run nothing on import.

Models (kind):
  'BO'      g channel, lam dU_g only (identical to the symmetric family: lam_c does not depend on r)
  'BOad'    g channel + K_gg - K_gg(inf) = dD_g + beta^2 (t_g - 1), counted from its own asymptote
  'BOad_thr' 'BOad' measured from the exact threshold (g asymptote lies +beta above it: no g/u mixing)
  'BOad+m'  'BOad' + the g-row mass functions (m_gg par/perp incl. beta^2 pieces), own asymptote
  '2ch'     channels (g, u): lam diag(dU_g, dU_u) + K' - 2 a d/dr - a', exact threshold (-beta)
  '2ch+m'   '2ch' + full 2x2 mass matrices (non-perturbative, Q = (1 + m/lam)^-1)   <- parameter-free theory
  '2ch+m_asym'  as '2ch+m' but m_uu, m_gu frozen at their atomic values (1/2 + 2b^2, -2b): error-bar variant
  '2ch+m_rc3'   as '2ch+m' with the u-row masses frozen below r = 3 instead of 2: error-bar variant

Numerics: regular solutions integrated outward from r0 = 1e-3 to r_x, two inward solutions (closed channel
decaying, open channel = zero-energy tail continued from RM = 2e4) from r_m = 40 to r_x; QR renormalization
in segments (sign of the determinant tracked); F(lam) = det[outward | inward] at r_x.  For one channel the
same code runs with a 1x1 system.  r = 1 (beta = 0) reproduces paper_data/bo_offset/coupling.py (tested).
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "solvers"))
sys.path.insert(0, os.path.join(ROOT, "paper_data", "bo_offset"))
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import coupling as C  # noqa: E402  (symmetric family: dU, dD, RB, RM, DHI)
from bo_h2plus import electronic as el  # noqa: E402
from bo_h2plus import unequal as UQ  # noqa: E402

R0, RX, RMATCH, RM = 1e-3, 6.0, 40.0, C.RM
RT = 60.0   # end of the g/u tables
RCUT_U = 2.0  # u-row masses frozen below this r (see Tables.masses)
MASS_KEYS = ["mv_gg0", "mv_gg2", "mv_uu0", "mv_uu2", "mv_gu1", "mr_gg0", "mr_gg2", "mr_uu0", "mr_uu2", "mr_gu1"]
MASS_ASYM = dict(mv_gg0=0.5, mv_gg2=2.0, mv_uu0=0.5, mv_uu2=2.0, mv_gu1=-2.0,
                 mr_gg0=0.5, mr_gg2=2.0, mr_uu0=0.5, mr_uu2=2.0, mr_gu1=-2.0)


# ---------------------------------------------------------------- tables
class Tables:
    """Splines of all beta-independent functions (unit mass, scaled r)."""

    def __init__(self, mass_csv=None, gu_path=None):
        tab = UQ.gu_table(path=gu_path) if gu_path else UQ.gu_table()
        R = tab["R"]
        x = np.log(R)
        self.Rmin, self.Rmax = R[0], R[-1]
        # r^2-weighted splines in log r (D_u ~ 2/r^2, k_gu ~ 1/r at small r)
        f = np.vstack([R * (tab["U_u"] + 0.5), R ** 2 * (tab["D_u"] - 0.25), tab["t_g_vir"] - 1.0,
                       tab["t_u_vir"] - 1.0, R * (tab["k_gu"] + 1.0), tab["dz_gu"]]).T
        self._gu = CubicSpline(x, f)
        self._gu_end = f[-1]
        self.mass = None
        path = mass_csv or os.path.join(HERE, "massfunc_gu.csv")
        if os.path.exists(path):
            import pandas as pd
            df = pd.read_csv(path).drop_duplicates("R").sort_values("R")
            self.mR = df.R.values
            self._m = CubicSpline(np.log(self.mR), df[MASS_KEYS].values)
            self.mass = df

    def gu(self, r):
        """(dU_u, dD_u, t_g-1, t_u-1, k_gu, dz_gu, d(dz_gu)/dr) at scalar r."""
        if r <= self.Rmin:
            v = self._gu(np.log(self.Rmin))
            rr = self.Rmin
            dU_u = (1.0 + (v[0] - 1.0) * r / rr) / r          # r dU_u -> 1 at r -> 0
            return dU_u, v[1] / r ** 2, v[2], v[3], v[4] / r - 1.0, v[5], 0.0
        if r >= self.Rmax:
            v = self._gu_end
            s4 = (self.Rmax / r) ** 4
            return (-9 / (4 * r ** 4) - 15 / (2 * r ** 6), v[1] / self.Rmax ** 2 * s4, -13.5 / r ** 4, -13.5 / r ** 4,
                    -1.0 + v[4] / self.Rmax * s4, 0.0, 0.0)
        lx = np.log(r)
        v = self._gu(lx)
        dv = self._gu(lx, 1) / r
        return v[0] / r, v[1] / r ** 2, v[2], v[3], v[4] / r - 1.0, v[5], dv[5]

    def masses(self, r, asym_ud=False, rcut_u=None):
        """dict of the mass-matrix coefficients at r (clamped below the grid, atomic values above).

        The u-row pieces (uu, gu) diverge as r -> 0 (2p sigma_u, 2p pi_u and 2s sigma_g become degenerate in
        the united atom); they are frozen at r < rcut_u (default RCUT_U) and multiplied by (r/rcut_u)^2 there (the u channel is
        strongly closed at such r, so any smooth choice is equivalent; variant '2ch+m_rc3' tests it)."""
        if self.mass is None:
            raise RuntimeError("massfunc_gu.csv not built (run build_tables.py)")
        if r >= self.mR[-1]:
            return dict(MASS_ASYM)
        v = self._m(np.log(max(r, self.mR[0])))
        d = dict(zip(MASS_KEYS, v))
        rc = RCUT_U if rcut_u is None else rcut_u
        if r < rc:
            # frozen at rc and switched off as (r/rc)^2 (keeps the r -> 0 Frobenius structure diagonal)
            vu = dict(zip(MASS_KEYS, self._m(np.log(rc))))
            for k in ("mv_uu0", "mv_uu2", "mv_gu1", "mr_uu0", "mr_uu2", "mr_gu1"):
                d[k] = vu[k] * (r / rc) ** 2
        if asym_ud:
            for k in ("mv_uu0", "mv_uu2", "mv_gu1", "mr_uu0", "mr_uu2", "mr_gu1"):
                d[k] = MASS_ASYM[k]
        return d


_TAB = None


def tables():
    global _TAB
    if _TAB is None:
        _TAB = Tables()
    return _TAB


# ---------------------------------------------------------------- model
class Model:
    def __init__(self, kind, r, L, tab=None):
        self.kind, self.r, self.L = kind, float(r), int(L)
        self.b = UQ.beta_of(self.r)
        self.tab = tab or tables()
        self.n = 2 if kind.startswith("2ch") else 1
        self.ad = kind != "BO"
        self.thr = kind == "BOad_thr"     # g channel only, measured from the exact threshold (asymptote +beta)
        self.withm = kind.endswith("+m") or kind.endswith("+m_asym")
        self.asym_ud = kind.endswith("_asym")
        self.rcut_u = 3.0 if kind.endswith("+m_rc3") else None
        if kind.endswith("+m_rc3"):
            self.withm = True
        if self.n == 2 and self.b == 0.0:
            raise ValueError("2-channel model needs r != 1 (at r = 1 g and u decouple; use 'BOad')")
        # matching radius: the closed (A-localized) channel must have decayed, kappa ~ sqrt(2 beta)
        self.rm = RMATCH if self.n == 1 else max(RMATCH, 30.0 / np.sqrt(2.0 * self.b))

    # -- coefficient matrices at r for coupling lam
    def coeffs(self, r, lam):
        """Return (Q_par^-1, Wmat, a, a') with  u' = Qinv w,  w' = -2 a u' - a' u + W u."""
        b, L = self.b, self.L
        dUg = float(C.dU(r))
        Kgg = float(C.dD(r)) if self.ad else 0.0
        dUu, dDu, tg1, tu1, kgu, dz, ddz = self.tab.gu(r)
        if self.ad:
            Kgg += b * b * tg1
        if self.thr:
            Kgg += b
        if self.n == 1:
            mpar = mperp = 0.0
            if self.withm:
                m = self.tab.masses(r)
                mpar, mperp = m["mv_gg0"] + b * b * m["mv_gg2"], m["mr_gg0"] + b * b * m["mr_gg2"]
            qinv = 1.0 + mpar / lam
            qperp = 1.0 / (1.0 + mperp / lam)
            W = np.array([[L * (L + 1) * qperp / r ** 2 + lam * dUg + Kgg]])
            return np.array([[qinv]]), W, None, None
        K = np.array([[Kgg + b, b * kgu], [b * kgu, dDu + b * b * tu1 + b]])
        a = np.array([[0.0, b * dz], [-b * dz, 0.0]])
        ap = np.array([[0.0, b * ddz], [-b * ddz, 0.0]])
        if self.withm:
            m = self.tab.masses(r, self.asym_ud, self.rcut_u)
            mpar = np.array([[m["mv_gg0"] + b * b * m["mv_gg2"], b * m["mv_gu1"]],
                             [b * m["mv_gu1"], m["mv_uu0"] + b * b * m["mv_uu2"]]])
            mperp = np.array([[m["mr_gg0"] + b * b * m["mr_gg2"], b * m["mr_gu1"]],
                              [b * m["mr_gu1"], m["mr_uu0"] + b * b * m["mr_uu2"]]])
            qinv = np.eye(2) + mpar / lam
            qperp = np.linalg.inv(np.eye(2) + mperp / lam)
        else:
            qinv, qperp = np.eye(2), np.eye(2)
        W = L * (L + 1) * qperp / r ** 2 + lam * np.diag([dUg, dUu]) + K
        return qinv, W, a, ap

    # -- open channel beyond r_m
    def open_channel(self, r, lam):
        """(q_par, q_perp, V) of the open channel at r >= RMATCH (lowest eigenvector of the potential matrix)."""
        L, b = self.L, self.b
        if self.n == 1:
            m0 = 0.5 + 2 * b * b if self.withm else 0.0
            q = lam / (lam + m0)
            if r <= C.RB:
                V = lam * float(C.dU(r)) + ((float(C.dD(r)) + b * b * self.tab.gu(r)[2]) if self.ad else 0.0)
            else:
                V = lam * float(C.dU(r)) + ((float(C.dD(r)) - b * b * 13.5 / r ** 4) if self.ad else 0.0)
            return q, q, V + L * (L + 1) * q / r ** 2
        m0 = 2.0 / (1.0 + self.r) ** 2 if self.withm else 0.0
        q = lam / (lam + m0)
        if r <= RT:
            _, W, _, _ = self.coeffs(r, lam)
            Wc = W - L * (L + 1) * self._qperp_cache(r, lam) / r ** 2
            V = float(np.linalg.eigvalsh(0.5 * (Wc + Wc.T))[0])
        else:
            # analytic tail: polarization + K' r^-4 (from the value at RT)
            _, W60, _, _ = self.coeffs(RT, lam)
            W60 = W60 - L * (L + 1) * self._qperp_cache(RT, lam) / RT ** 2
            v60 = float(np.linalg.eigvalsh(0.5 * (W60 + W60.T))[0]) - lam * float(C.dU(RT))
            V = lam * float(C.dU(r)) + v60 * (RT / r) ** 4
        return q, q, V + L * (L + 1) * q / r ** 2

    def _qperp_cache(self, r, lam):
        if not self.withm:
            return np.eye(2)
        m = self.tab.masses(r, self.asym_ud, self.rcut_u)
        b = self.b
        mperp = np.array([[m["mr_gg0"] + b * b * m["mr_gg2"], b * m["mr_gu1"]],
                          [b * m["mr_gu1"], m["mr_uu0"] + b * b * m["mr_uu2"]]])
        return np.linalg.inv(np.eye(2) + mperp / lam)

    # -- ODE right-hand side for a set of column solutions y = [u (n x k); w (n x k)]
    def rhs(self, lam, k):
        n = self.n

        def f(r, y):
            Y = y.reshape(2 * n, k)
            U, Wv = Y[:n], Y[n:]
            qinv, W, a, ap = self.coeffs(r, lam)
            Up = qinv @ Wv
            Wp = W @ U
            if a is not None:
                Wp = Wp - 2 * a @ Up - ap @ U
            return np.vstack([Up, Wp]).ravel()
        return f

    def tail(self, lam):
        """Open-channel zero-energy solution with no growing part, integrated inward from RM to RMATCH: (u, w)."""
        L = self.L
        q, _, _ = self.open_channel(RM, lam)
        Vrm = self.open_channel(RM, lam)[2] - L * (L + 1) * q / RM ** 2
        if L == 0:
            C4 = max(-Vrm * RM ** 4, 1e-30)
            bb = np.sqrt(C4 / q)
            u, up = RM * np.sin(bb / RM), np.sin(bb / RM) - (bb / RM) * np.cos(bb / RM)
        else:
            u, up = RM ** (-L), -L * RM ** (-L - 1)

        def f(r, y):
            q, _, V = self.open_channel(r, lam)
            return [y[1] / q, V * y[0]]
        s = solve_ivp(f, (RM, self.rm), [u, q * up], method="DOP853", rtol=1e-11, atol=1e-300)
        return s.y[0, -1], s.y[1, -1]

    def F(self, lam):
        """Matching determinant (sign-continuous, magnitude normalized)."""
        n, L = self.n, self.L
        # outward: regular solutions u ~ r^(L+1) (1 + lam r/(2L+2)) in each channel
        a1 = lam / (2 * L + 2)
        Y = np.zeros((2 * n, n))
        qinv0, _, _, _ = self.coeffs(R0, lam)
        for i in range(n):
            Y[i, i] = R0 ** (L + 1) * (1 + a1 * R0)
            up = (L + 1) * R0 ** L + (L + 2) * a1 * R0 ** (L + 1)
            Y[n:, i] = np.linalg.solve(qinv0, np.eye(n)[:, i] * up)
        sgn = 1.0
        Yo, s1 = self._propagate(lam, Y, np.geomspace(R0, RX, 25))
        sgn *= s1
        # inward: closed decaying + open tail at RMATCH
        RMm = self.rm
        qinv, W, _, _ = self.coeffs(RMm, lam)
        ut, wt = self.tail(lam) if not self.thr else (None, None)
        if n == 1 and self.thr:
            # g channel closed at the exact threshold: decaying solution exp(-kappa r), kappa^2 = V(r_m)/q
            q, _, V = self.open_channel(RMm, lam)
            V += self.b
            kap = np.sqrt(max(V, 1e-14) / q)
            Yi = np.array([[1.0], [-q * kap]])
        elif n == 1:
            Yi = np.array([[ut], [wt]])
        else:
            Wc = W - L * (L + 1) * self._qperp_cache(RMm, lam) / RMm ** 2
            ev, V = np.linalg.eigh(0.5 * (Wc + Wc.T))
            vo, vc = V[:, 0], V[:, 1]
            if vo.sum() < 0:
                vo = -vo
            if vc[0] < 0:
                vc = -vc
            qo = 1.0 / float(vo @ qinv @ vo)
            qc = 1.0 / float(vc @ qinv @ vc)
            kap = np.sqrt(max(ev[1] - L * (L + 1) * qc / RMm ** 2, 1e-12) / qc)
            # open: u = vo ut, u' = vo ut'/qo -> w = Q u'
            Q = np.linalg.inv(qinv)
            uo, uo_p = vo * ut, vo * (wt / qo)
            uc, uc_p = vc, -kap * vc
            Yi = np.column_stack([np.concatenate([uo, Q @ uo_p]), np.concatenate([uc, Q @ uc_p])])
        Yi, s2 = self._propagate(lam, Yi, np.geomspace(RMm, RX, 40))
        sgn *= s2
        M = np.hstack([Yo, Yi])
        M = M / np.linalg.norm(M, axis=0)
        return sgn * np.linalg.det(M)

    def _propagate(self, lam, Y, pts):
        n, k = self.n, Y.shape[1]
        f = self.rhs(lam, k)
        sgn = 1.0
        for a, b in zip(pts[:-1], pts[1:]):
            s = solve_ivp(f, (a, b), Y.ravel(), method="DOP853", rtol=1e-11, atol=1e-30)
            if not s.success:
                raise RuntimeError(s.message)
            Y = s.y[:, -1].reshape(2 * n, k)
            Qm, Rm = np.linalg.qr(Y)
            sgn *= np.sign(np.linalg.det(Rm))
            Y = Qm
        return Y, sgn

    def roots(self, lo, hi, n=40):
        lams = np.linspace(lo, hi, n)
        Fv = np.array([self.F(x) for x in lams])
        out = []
        for i in range(n - 1):
            if np.sign(Fv[i]) != np.sign(Fv[i + 1]):
                out.append(brentq(self.F, lams[i], lams[i + 1], xtol=1e-10, rtol=1e-13))
        return out, lams, Fv


def critical(kind, r, L, lo, hi, n=40, tab=None):
    """All lam_c of the model in [lo, hi] and the corresponding M_c (light-nucleus mass)."""
    m = Model(kind, r, L, tab)
    lams, _, _ = m.roots(lo, hi, n)
    return [(lam, UQ.M_of_lam(lam, r)) for lam in lams]
