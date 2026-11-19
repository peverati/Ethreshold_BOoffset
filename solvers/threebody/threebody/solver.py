"""Generalized eigenproblem, eigenvectors and observables."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import eigh
from scipy.special import roots_laguerre, roots_legendre

from . import laguerre
from .basis import PerimetricBasis
from .operators import PAIRS, kinetic_terms, overlap_terms, power_terms, third
from .system import ThreeBodySystem


class MatrixSet:
    """Mass- and charge-independent matrices of one basis: S, T_i (per unit 1/m_i), V_ij.

    H(system) = sum_i (1/m_i) T_i + sum_ij Z_ij V_ij, so a mass scan at fixed
    basis only re-diagonalizes.
    """

    def __init__(self, basis: PerimetricBasis):
        self.basis = basis
        cache = {}
        t0 = time.perf_counter()
        self.S = basis.matrix(overlap_terms(), cache)
        self.T = [basis.matrix(kinetic_terms(i), cache) for i in range(3)]
        self.V = {p: basis.matrix(power_terms(p[0], p[1], -1), cache) for p in PAIRS}
        self.build_time = time.perf_counter() - t0
        self._orth = None

    def hamiltonian(self, system: ThreeBodySystem):
        im = system.inverse_masses
        H = sum(im[i] * self.T[i] for i in range(3) if im[i] != 0.0)
        for p, z in system.coupling_dict().items():
            if z != 0.0:
                H = H + z * self.V[p]
        return H

    def orthogonalizer(self, tol=1e-13):
        """Canonical orthogonalization X (X^T S X = 1), dropping near-null directions of S.

        S is first scaled to unit diagonal; eigenvectors with eigenvalue below
        tol * (largest eigenvalue) are discarded.
        """
        if self._orth is None or self._orth[0] != tol:
            d = 1.0 / np.sqrt(np.diag(self.S))
            Ss = self.S * np.outer(d, d)
            lam, U = np.linalg.eigh(Ss)
            keep = lam > tol * lam[-1]
            X = (d[:, None] * U[:, keep]) / np.sqrt(lam[keep])
            self._orth = (tol, X, lam[0], lam[-1], int((~keep).sum()))
        return self._orth


@dataclass
class Result:
    system: ThreeBodySystem
    basis: PerimetricBasis
    energies: np.ndarray
    vectors: np.ndarray  # columns, normalized c^T S c = 1, reduced (symmetry-adapted) basis
    info: dict = field(default_factory=dict)
    matrices: MatrixSet = None

    @property
    def threshold(self):
        return self.system.threshold()

    def n_below_threshold(self, margin=1e-10):
        """Number of variational roots below E_th - margin.

        Each such root is a rigorous upper bound (Hylleraas-Undheim-MacDonald),
        so this is a *lower bound* on the number of bound states of this
        symmetry, up to rounding (controlled by margin).
        """
        return int(np.sum(self.energies < self.threshold - margin))

    # ------------------------------------------------------------- wave function
    def tensor(self, state=0):
        return self.basis.full_tensor(self.vectors[:, state])

    def expectation_terms(self, terms, state=0, other=None):
        """<psi|O|psi> for a term list, via exact tensor contraction (no N^2 matrices)."""
        C = self.tensor(state)
        Cb = C if other is None else other
        b = self.basis
        cache = {}
        total = 0.0
        for coef, p, k, l in terms:
            Ms = []
            for m in range(3):
                key = (m, p[m], int(k == m), int(l == m))
                if key not in cache:
                    cache[key] = laguerre.one_d_matrix(b.nmax[m], p[m], key[2], key[3], b.alphas[m])
                Ms.append(cache[key])
            X = np.einsum("abc,ad->dbc", C, Ms[0], optimize=True)
            X = np.einsum("dbc,be->dec", X, Ms[1], optimize=True)
            X = np.einsum("dec,cf->def", X, Ms[2], optimize=True)
            total += coef * np.sum(X * Cb)
        return total

    def norm(self, state=0):
        return self.expectation_terms(overlap_terms(), state)

    def expect_power(self, i, j, n, state=0):
        """<r_ij^n> for n >= -1 (exact)."""
        return self.expectation_terms(power_terms(i, j, n), state) / self.norm(state)

    def energy_components(self, state=0):
        im = self.system.inverse_masses
        nrm = self.norm(state)
        T = sum(im[i] * self.expectation_terms(kinetic_terms(i), state) for i in range(3) if im[i] != 0) / nrm
        V = 0.0
        for p, z in self.system.coupling_dict().items():
            if z != 0.0:
                V += z * self.expectation_terms(power_terms(p[0], p[1], -1), state) / nrm
        return dict(T=T, V=V, E=T + V, virial_ratio=-V / T)

    def dE_dinverse_mass(self, i, state=0):
        """Hellmann-Feynman derivative dE/d(1/m_i) = <T_i> at fixed basis (exact for the
        variational eigenvalue, because the basis does not depend on the masses).
        For the {M,M,1} family: dE/dM = -(<T_0> + <T_1>) / M^2."""
        return self.expectation_terms(kinetic_terms(i), state) / self.norm(state)

    # -------------------------------------------------- pair distribution function
    def pair_distribution(self, i, j, R, state=0, nq=None):
        """P_ij(R): probability density of r_ij (int P dR = 1).

        g(R) = P(R) / (4 pi R^2) satisfies int 4 pi R^2 g dR = 1.
        Fix r_ij = R, i.e. s_i + s_j = 2R; with s_i = t in [0, 2R], s_k = u >= 0,
        dr_ik dr_jk = dt du / 2 and the measure is R r_ik r_jk.
        The u integral is done exactly with 1D moment matrices, the t integral
        by Gauss-Legendre (exact for alpha_i = alpha_j, converged otherwise).
        """
        i, j = min(i, j), max(i, j)
        k = third(i, j)
        b = self.basis
        C = self.tensor(state)
        R = np.atleast_1d(np.asarray(R, dtype=float))
        if nq is None:
            nq = b.nmax[i] + b.nmax[j] + 30
        xg, wg = roots_legendre(nq)
        Sk = [laguerre.moment_matrix(b.nmax[k], p, b.alphas[k]) for p in range(3)]
        # move axes so that C has order (i, j, k)
        Cijk = np.transpose(C, (i, j, k))
        out = np.empty_like(R)
        nrm = self.norm(state)
        for n, RR in enumerate(R):
            if RR <= 0:
                out[n] = 0.0
                continue
            t = RR * (xg + 1.0)
            w = RR * wg
            Fi = laguerre.evaluate(b.nmax[i], b.alphas[i], t)  # (ni, q)
            Fj = laguerre.evaluate(b.nmax[j], b.alphas[j], 2 * RR - t)
            # v[q, c] = sum_ab C[a,b,c] Fi[a,q] Fj[b,q]
            v = np.einsum("abc,aq,bq->qc", Cijk, Fi, Fj, optimize=True)
            # u-integral of |psi|^2 (t(2R-t) + 2R u + u^2)/4
            q0 = np.einsum("qc,cd,qd->q", v, Sk[0], v)
            q1 = np.einsum("qc,cd,qd->q", v, Sk[1], v)
            q2 = np.einsum("qc,cd,qd->q", v, Sk[2], v)
            integrand = (t * (2 * RR - t) * q0 + 2 * RR * q1 + q2) / 4.0
            out[n] = RR * 0.5 * np.sum(w * integrand)
        return out / nrm

    def g_of_r(self, i, j, R, state=0):
        R = np.asarray(R, dtype=float)
        P = self.pair_distribution(i, j, R, state)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.where(R > 0, P / (4 * np.pi * R**2), np.nan)
        return g

    def radial_moment_quadrature(self, i, j, n, state=0, nq=None):
        """<r_ij^n> from Gauss-Laguerre quadrature of P_ij(R) (works for n = -2)."""
        b = self.basis
        i, j = min(i, j), max(i, j)
        gam = b.alphas[i] + b.alphas[j]
        if nq is None:
            nq = b.nmax[i] + b.nmax[j] + 20
        x, w = roots_laguerre(nq)
        R = x / gam
        P = self.pair_distribution(i, j, R, state)
        return float(np.sum(w * np.exp(x) / gam * P * R**n))

    # ------------------------------------------------------------- observables
    def observables(self, state=0, inverse_square=True):
        obs = {}
        for (i, j) in PAIRS:
            tag = f"{i}{j}"
            obs[f"r{tag}"] = self.expect_power(i, j, 1, state)
            obs[f"r{tag}^2"] = self.expect_power(i, j, 2, state)
            obs[f"1/r{tag}"] = self.expect_power(i, j, -1, state)
            if inverse_square:
                obs[f"1/r{tag}^2"] = self.radial_moment_quadrature(i, j, -2, state)
        ec = self.energy_components(state)
        obs.update(ec)
        r, r2 = obs["r01"], obs["r01^2"]
        obs["lindemann01"] = np.sqrt(max(r2 - r * r, 0.0)) / r
        return obs

    # ---------------------------------------------------------------- persistence
    def save(self, path, states=None):
        """Save basis definition, system, energies and coefficients (npz + JSON header)."""
        if states is None:
            states = range(len(self.energies))
        states = list(states)
        meta = dict(system=dict(masses=list(self.system.masses), charges=list(self.system.charges),
                                couplings={f"{k[0]}{k[1]}": v for k, v in self.system.coupling_dict().items()},
                                identical=self.system.identical, name=self.system.name),
                    basis=self.basis.describe(), info=_jsonable(self.info),
                    convention="perimetric s_k = r_kj + r_kl - r_jl; f_n(s)=L_n(alpha s)exp(-alpha s/2); "
                               "psi = sum_a C[a0,a1,a2] f_a0(s0) f_a1(s1) f_a2(s2); measure r01 r02 r12 ds/4")
        np.savez(path, meta=json.dumps(meta), energies=self.energies[states],
                 vectors=self.vectors[:, states], index=self.basis.index,
                 tensors=np.stack([self.tensor(s) for s in states]))


def _jsonable(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        elif isinstance(v, np.ndarray):
            v = v.tolist()
        out[k] = v
    return out


def solve(system: ThreeBodySystem, basis: PerimetricBasis, nstates=5, matrices: MatrixSet = None,
          tol=1e-13) -> Result:
    """Lowest nstates eigenpairs of the given symmetry."""
    if system.identical is None and basis.symmetry is not None:
        raise ValueError("symmetry-adapted basis requires identical particles 0 and 1")
    if matrices is None:
        matrices = MatrixSet(basis)
    t0 = time.perf_counter()
    H = matrices.hamiltonian(system)
    _, X, lmin, lmax, ndrop = matrices.orthogonalizer(tol)
    Hp = X.T @ H @ X
    k = min(nstates, Hp.shape[0])
    e, U = eigh(Hp, subset_by_index=[0, k - 1])
    C = X @ U
    # deterministic sign: <F_000|psi> > 0 if not tiny, else largest component positive
    SC = matrices.S @ C
    for n in range(k):
        ref = SC[np.argmax(np.abs(SC[:, n])), n] if abs(SC[0, n]) < 1e-6 else SC[0, n]
        if ref < 0:
            C[:, n] *= -1
    info = dict(size=basis.size, kept=X.shape[1], dropped=ndrop, S_min_eig_scaled=lmin,
                S_cond_scaled=lmax / lmin if lmin > 0 else np.inf,
                build_time=matrices.build_time, solve_time=time.perf_counter() - t0)
    return Result(system, basis, e, C, info, matrices)


def overlap_between(res_a: Result, res_b: Result, state_a=0, state_b=0, normalized=True):
    """<psi_a|psi_b> in fixed physical (interparticle) coordinates.

    The two results may use different scale parameters alpha and sizes;
    the 1D integrals between different scales are evaluated exactly by
    Gauss-Laguerre quadrature.  Both must describe the same particle labelling.
    """
    Ca, Cb = res_a.tensor(state_a), res_b.tensor(state_b)
    ba, bb = res_a.basis, res_b.basis
    cache = {}
    total = 0.0
    for coef, p, k, l in overlap_terms():
        Ms = []
        for m in range(3):
            key = (m, p[m])
            if key not in cache:
                cache[key] = laguerre.mixed_moment_matrix(ba.nmax[m], ba.alphas[m], bb.nmax[m], bb.alphas[m], p[m])
            Ms.append(cache[key])
        X = np.einsum("abc,ad,be,cf->def", Ca, Ms[0], Ms[1], Ms[2], optimize=True)
        total += coef * np.sum(X * Cb)
    if normalized:
        total /= np.sqrt(res_a.norm(state_a) * res_b.norm(state_b))
    return total
