"""Double-double (~32 digit) assembly of the coupled-basis matrices (high-L precision fix).  Written 2026-11-17.

Why: at L = 7 the angular term lists have up to ~1000 terms with powers up to s^16; each matrix element is a sum
of large monomial integrals that cancel, so double-precision assembly (threebody.coupled.CoupledBasis._assemble)
leaves an absolute error far above the small eigenvalues of S (validate_L 7 fails; spurious roots at tol 1e-14).
Here the same term lists are evaluated with
  * exact 1D factors: same-scale integrals int f_a f_b s^p ds = alpha^(da+db-p-1) * Q / 2^(da+db) with Q an exact
    integer matrix ((2D)^T J^p (2D) in Python integers, J = Laguerre Jacobi matrix, D = derivative matrix);
  * double-double accumulation of coef * X0 * X1 * X2 over all terms (error-free transformations, vectorized).
The result (hi + lo) is rounded to double; projection and eigensolution stay in the existing double code.
Mixed-scale blocks (alpha_a != alpha_b in some coordinate) use mixed_moment_dd (exact rational expansion evaluated
with mpmath at MP_DPS digits, then split to double-double).

Nothing in solvers/ is modified: HPCoupledBasis subclasses CoupledBasis and only overrides _assemble.
"""
from __future__ import annotations

import os
from fractions import Fraction

import numpy as np

import tb_hp_local  # noqa: F401  (private term cache; must precede threebody imports)
from threebody import angular
from threebody.coupled import CoupledBasis

MP_DPS = 80
_SPLIT = 134217729.0  # 2^27 + 1 (Dekker)


# ------------------------------------------------------------------ double-double primitives (vectorized)
def two_sum(a, b):
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)


def quick_two_sum(a, b):
    s = a + b
    return s, b - (s - a)


def split(a):
    t = _SPLIT * a
    hi = t - (t - a)
    return hi, a - hi


def two_prod(a, b):
    p = a * b
    ah, al = split(a)
    bh, bl = split(b)
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


def dd_add(ah, al, bh, bl):
    s, e = two_sum(ah, bh)
    e = e + (al + bl)
    return quick_two_sum(s, e)


def dd_mul(ah, al, bh, bl):
    p, e = two_prod(ah, bh)
    e = e + (ah * bl + al * bh)
    return quick_two_sum(p, e)


def to_dd(x):
    """Exact Python int / Fraction / mpf (array of objects or scalar) -> (hi, lo) float64 arrays."""
    arr = np.asarray(x, dtype=object)
    hi = np.vectorize(float, otypes=[float])(arr) if arr.size else np.zeros(arr.shape)
    lo = np.vectorize(lambda v, h: float(Fraction(v) - Fraction(h)) if not hasattr(v, "_mpf_") else float(v - h),
                      otypes=[float])(arr, hi) if arr.size else np.zeros(arr.shape)
    return hi, lo


def dd_scalar_pow(a, n):
    """a^n (a double, n integer, possibly negative) as an exact-to-dd (hi, lo) pair via exact Fractions."""
    v = Fraction(a) ** n
    hi = float(v)
    return hi, float(v - Fraction(hi))


# ------------------------------------------------------------------ exact 1D factors
_qcache = {}


def _q_same(nmax_a, nmax_b, p, da, db):
    """Exact integer matrix Q with int (d^da f_a)(d^db f_b) s^p ds = alpha^(da+db-p-1) Q / 2^(da+db)."""
    key = (nmax_a, nmax_b, p, da, db)
    if key in _qcache:
        return _qcache[key]
    n = max(nmax_a, nmax_b)
    size = n + 1 + p
    J = np.zeros((size, size), dtype=object)
    for i in range(size):
        J[i, i] = 2 * i + 1
        if i + 1 < size:
            J[i, i + 1] = J[i + 1, i] = -(i + 1)
    Jp = np.eye(size, dtype=int).astype(object)
    for _ in range(p):
        Jp = Jp.dot(J)
    Q = Jp[: n + 1, : n + 1]
    D2 = np.zeros((n + 1, n + 1), dtype=object)       # 2 D: -2 above the diagonal, -1 on it
    for i in range(n + 1):
        D2[i, i] = -1
        for j in range(i + 1, n + 1):
            D2[i, j] = -2
    if da:
        Q = D2.T.dot(Q)
    if db:
        Q = Q.dot(D2)
    Q = Q[: nmax_a + 1, : nmax_b + 1]
    _qcache[key] = Q
    return Q


_glcache = {}


def gauss_laguerre_mp(nq):
    """Gauss-Laguerre nodes/weights (weight exp(-x)) to MP_DPS digits: scipy double roots refined by Newton."""
    import mpmath as mp
    from scipy.special import roots_laguerre
    if nq in _glcache:
        return _glcache[nq]
    mp.mp.dps = MP_DPS

    def lag(n, x):          # L_n(x), L_{n-1}(x) (standard, L_n(0) = 1)
        l0, l1 = mp.mpf(1), 1 - x
        if n == 0:
            return l0, mp.mpf(0)
        for k in range(1, n):
            l0, l1 = l1, ((2 * k + 1 - x) * l1 - k * l0) / (k + 1)
        return l1, l0

    xs, ws = [], []
    for x0 in roots_laguerre(nq)[0]:
        x = mp.mpf(x0)
        for _ in range(100):
            ln, lm = lag(nq, x)
            dx = ln / (nq * (ln - lm) / x)      # L_n' = n (L_n - L_{n-1}) / x
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-MP_DPS + 5) * (1 + abs(x)):
                break
        l1, _ = lag(nq + 1, x)
        xs.append(x)
        ws.append(x / ((nq + 1) ** 2 * l1 ** 2))
    _glcache[nq] = (xs, ws)
    return xs, ws


def mixed_moment_dd(nmax_a, alpha, nmax_b, beta, p, da, db):
    """int (d^da f_a(s;alpha)) (d^db f_b(s;beta)) s^p ds, alpha != beta: Gauss-Laguerre in mpmath (exact for the
    polynomial x exp(-(alpha+beta)s/2) integrand up to MP_DPS rounding), derivatives via exact D."""
    import mpmath as mp
    mp.mp.dps = MP_DPS
    a = mp.mpf(Fraction(alpha).numerator) / Fraction(alpha).denominator
    b = mp.mpf(Fraction(beta).numerator) / Fraction(beta).denominator
    g = (a + b) / 2
    nq = (nmax_a + nmax_b + p) // 2 + 2
    xs, ws = gauss_laguerre_mp(nq)

    def fvals(nmax, al, s):     # f_n(s; al) without the exponential (it is in the weight): L_n(al s)
        x = al * s
        out = [mp.mpf(1), 1 - x]
        for n in range(1, nmax):
            out.append(((2 * n + 1 - x) * out[n] - n * out[n - 1]) / (n + 1))
        return out[: nmax + 1]

    # exp(-a s/2) exp(-b s/2) = exp(-g s) = exp(-x): the Gauss-Laguerre weight; s = x / g, ds = dx / g
    Fa = np.empty((nmax_a + 1, nq), dtype=object)
    Fb = np.empty((nmax_b + 1, nq), dtype=object)
    for q, (x, w) in enumerate(zip(xs, ws)):
        sq = x / g
        Wq = w / g * sq ** p
        fa, fb = fvals(nmax_a, a, sq), fvals(nmax_b, b, sq)
        for i in range(nmax_a + 1):
            Fa[i, q] = fa[i] * Wq
        for j in range(nmax_b + 1):
            Fb[j, q] = fb[j]
    S = Fa.dot(Fb.T)
    if da:   # alpha * (D^T S), D_kn = -1 (k < n), -1/2 (k = n)
        cum = np.cumsum(np.vstack([np.zeros((1, S.shape[1]), dtype=object), S[:-1]]), axis=0)
        S = a * (-S / 2 - cum)
    if db:
        cum = np.cumsum(np.hstack([np.zeros((S.shape[0], 1), dtype=object), S[:, :-1]]), axis=1)
        S = b * (-S / 2 - cum)
    return to_dd(S)


def mixed_moment_exact_mp(nmax_a, alpha, nmax_b, beta, p, da, db):
    """Reference (slow, small n only; used in the tests): explicit power expansion of the Laguerre polynomials.

    int (d^da f_a(s;alpha)) (d^db f_b(s;beta)) s^p ds for alpha != beta, to ~1e-30 relative (mpmath).

    f_n(s; a) = L_n(a s) exp(-a s/2) with L_n(x) = sum_k (-1)^k C(n,k) x^k / k!; derivatives via
    d f_n/ds = a sum_k D_kn f_k.  int s^m exp(-g s) ds = m! / g^(m+1), g = (alpha + beta)/2."""
    import mpmath as mp
    mp.mp.dps = MP_DPS
    a, b = mp.mpf(Fraction(alpha).numerator) / Fraction(alpha).denominator, \
        mp.mpf(Fraction(beta).numerator) / Fraction(beta).denominator
    g = (a + b) / 2
    na, nb = nmax_a, nmax_b
    m_max = na + nb + p
    mom = [mp.factorial(m) / g ** (m + 1) for m in range(m_max + 1)]
    ca = [[(-1) ** k * mp.binomial(n, k) / mp.factorial(k) * a ** k for k in range(n + 1)] for n in range(na + 1)]
    cb = [[(-1) ** k * mp.binomial(n, k) / mp.factorial(k) * b ** k for k in range(n + 1)] for n in range(nb + 1)]
    S = np.empty((na + 1, nb + 1), dtype=object)
    for i in range(na + 1):
        for j in range(nb + 1):
            S[i, j] = mp.fsum(ca[i][k] * cb[j][l] * mom[p + k + l] for k in range(i + 1) for l in range(j + 1))
    if da:   # alpha * (D^T S), D_kn = -1 (k < n), -1/2 (k = n)
        T = np.empty_like(S)
        for i in range(na + 1):
            for j in range(nb + 1):
                T[i, j] = a * (-mp.mpf(1) / 2 * S[i, j] - mp.fsum(S[k, j] for k in range(i)))
        S = T
    if db:
        T = np.empty_like(S)
        for i in range(na + 1):
            for j in range(nb + 1):
                T[i, j] = b * (-mp.mpf(1) / 2 * S[i, j] - mp.fsum(S[i, k] for k in range(j)))
        S = T
    return to_dd(S)


_odcache = {}


def one_d_dd(alpha_a, na, alpha_b, nb, p, da, db):
    """(hi, lo) of the 1D matrix int (d^da f_a)(d^db f_b) s^p ds, a <= na, b <= nb."""
    key = (alpha_a, na, alpha_b, nb, p, da, db)
    if key in _odcache:
        return _odcache[key]
    if alpha_a == alpha_b:
        Q = _q_same(na, nb, p, da, db)
        qh, ql = to_dd(Q)
        sh, sl = dd_scalar_pow(alpha_a, da + db - p - 1)
        sh, sl = sh / 2 ** (da + db), sl / 2 ** (da + db)        # exact (power of two)
        out = dd_mul(qh, ql, sh, sl)
    else:
        out = mixed_moment_dd(na, alpha_a, nb, alpha_b, p, da, db)
    _odcache[key] = out
    return out


# ------------------------------------------------------------------ raw block and assembly
def raw_matrix_dd(terms, A, B):
    """Double-double version of threebody.threshold.raw_matrix: returns (hi, lo)."""
    na, nb = A.nmax, B.nmax
    od = lambda m, p, da, db: one_d_dd(A.alphas[m], na[m], B.alphas[m], nb[m], p, da, db)
    hi = np.zeros((A.size, B.size))
    lo = np.zeros((A.size, B.size))
    r = [A.index[:, m][:, None] for m in range(3)]
    c = [B.index[:, m][None, :] for m in range(3)]
    for coef, p, k, l in terms:
        f = [(p[m], int(k == m), int(l == m)) for m in range(3)]
        x0h, x0l = od(0, *f[0])
        x1h, x1l = od(1, *f[1])
        x2h, x2l = od(2, *f[2])
        th, tl = dd_mul(x0h[r[0], c[0]], x0l[r[0], c[0]], x1h[r[1], c[1]], x1l[r[1], c[1]])
        th, tl = dd_mul(th, tl, x2h[r[2], c[2]], x2l[r[2], c[2]])
        th, tl = dd_mul(th, tl, float(coef), 0.0)
        hi, lo = dd_add(hi, lo, th, tl)
    return hi, lo


# ------------------------------------------------------------------ parallel assembly (fork workers)
NPROC = int(os.environ.get("HP_NPROC", "20"))
_CTX = {}


def _od_keys(terms, A, B):
    na, nb = A.nmax, B.nmax
    out = set()
    for coef, p, k, l in terms:
        for m in range(3):
            out.add((A.alphas[m], na[m], B.alphas[m], nb[m], p[m], int(k == m), int(l == m)))
    return out


def _od_worker(key):
    return key, one_d_dd(*key)


def _block_worker(IJ):
    I, J = IJ
    units, op, eps, sym = _CTX["units"], _CTX["op"], _CTX["eps"], _CTX["sym"]
    c, B, nb = units[I]
    d, C, nc = units[J]
    bh, bl = raw_matrix_dd(angular.operator_terms(op, c, d), B, C)
    if sym is not None:
        mh, ml = raw_matrix_dd(angular.operator_terms(op, c, angular.swap_component(d)), B, C.mirror())
        bh, bl = dd_add(bh, bl, eps * mh, eps * ml)
        w = 2.0 * np.outer(nb, nc)
        bh, bl = dd_mul(bh, bl, w, np.zeros_like(w))
    return I, J, bh, bl


def _pairs(units):
    return [(I, J) for I in range(len(units)) for J in range(I, len(units))]


class HPCoupledBasis(CoupledBasis):
    """CoupledBasis with double-double assembly (identical interface; matrices rounded to double at the end).
    Assembly runs on HP_NPROC fork workers (1D factors precomputed in parallel first, inherited by the block
    workers).  self.mats_lo keeps the low parts (diagnostics)."""

    two_pass = os.environ.get("HP_TWO_PASS", "0") == "1"

    def project(self, tol=1e-12):
        """NOT USED FOR RESULTS (2026-11-17 09:15): the two-pass option did not pass validate_hp (L = 7 tol 1e-14
        non-monotone by 3.8e-8 at N = 20; tol 1e-15 collapses for N <= 16), because X^T S X itself is rounding-limited
        in double at these tols (measured ortho_err ~1e-4 = eps x condition).  Kept for the record; default off.
        Canonical orthogonalization; with two_pass (HP_TWO_PASS=1) a second pass re-orthogonalizes X^T S X and
        re-applies the threshold (single-pass X is not orthonormal for tol <= 1e-14 at L = 7: X^T S X eigenvalues
        0.70..1.22 at tol 1e-14, 0.04..2.4 at 1e-15, which produced spurious levels; diag_spurious.txt).  Records
        ortho_err = max|X^T S X - I| of the final X."""
        if self._proj is not None and self._proj["tol"] == tol:
            return self._proj
        pr = CoupledBasis.project(self, tol)
        S = self.mats["S"]
        if self.two_pass:
            X = pr["X"]
            lam2, U2 = np.linalg.eigh(X.T @ S @ X)
            keep2 = lam2 > tol * lam2[-1]
            X = X @ (U2[:, keep2] / np.sqrt(lam2[keep2]))
            pr.update(X=X, kept=int(keep2.sum()), dropped=pr["dropped"] + int((~keep2).sum()),
                      pass2_min=float(lam2[0]), pass2_max=float(lam2[-1]), P={})
        pr["ortho_err"] = float(np.max(np.abs(pr["X"].T @ S @ pr["X"] - np.eye(pr["X"].shape[1]))))
        self._proj = pr
        return pr

    def _assemble(self, op):
        import multiprocessing as mpr
        units = self.units
        pairs = _pairs(units)
        keys = set()
        for I, J in pairs:
            c, B, _ = units[I]
            d, C, _ = units[J]
            keys |= _od_keys(angular.operator_terms(op, c, d), B, C)
            if self.symmetry is not None:
                keys |= _od_keys(angular.operator_terms(op, c, angular.swap_component(d)), B, C.mirror())
        todo = sorted(k for k in keys if k not in _odcache)
        ctx = mpr.get_context("fork")
        if todo:
            if NPROC > 1 and len(todo) > 1:
                with ctx.Pool(min(NPROC, len(todo))) as pool:
                    for k, v in pool.imap_unordered(_od_worker, todo):
                        _odcache[k] = v
            else:
                for k in todo:
                    _odcache[k] = one_d_dd(*k)
        n = self.size
        hi = np.zeros((n, n))
        lo = np.zeros((n, n))
        off = self.offsets
        _CTX.update(units=units, op=op, eps=self.eps, sym=self.symmetry)

        def put(I, J, bh, bl):
            hi[off[I]:off[I + 1], off[J]:off[J + 1]] = bh
            lo[off[I]:off[I + 1], off[J]:off[J + 1]] = bl
            if J != I:
                hi[off[J]:off[J + 1], off[I]:off[I + 1]] = bh.T
                lo[off[J]:off[J + 1], off[I]:off[I + 1]] = bl.T

        if NPROC > 1 and len(pairs) > 1:
            with ctx.Pool(min(NPROC, len(pairs))) as pool:
                for res in pool.imap_unordered(_block_worker, pairs):
                    put(*res)
        else:
            for IJ in pairs:
                put(*_block_worker(IJ))
        if not hasattr(self, "mats_lo"):
            self.mats_lo = {}
        sh, sl = dd_add(hi, lo, hi.T, lo.T)
        self.mats_lo[op] = 0.5 * sl
        return 0.5 * (sh + sl)
