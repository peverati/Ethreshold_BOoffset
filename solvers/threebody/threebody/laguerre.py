"""One-dimensional Laguerre machinery for the perimetric basis.

The 1D basis functions are

    f_n(s; alpha) = L_n(alpha s) exp(-alpha s / 2),   n = 0, 1, 2, ...

with L_n the (orthonormal w.r.t. exp(-x)) Laguerre polynomials.  All
one-dimensional integrals needed by the solver are of the form

    M[a, b] = int_0^inf (d^da f_a)(s) (d^db f_b)(s) s^p ds ,   da, db in {0, 1}

and are obtained *exactly* (integer arithmetic up to powers of alpha) from

* the Jacobi matrix J of the Laguerre polynomials (x L_n = -(n+1) L_{n+1}
  + (2n+1) L_n - n L_{n-1}), so that int f_a f_b s^p ds = alpha^-(p+1) (J^p)_ab;
* the derivative matrix D, d/ds f_n = alpha sum_k D_kn f_k with
  D_kn = -1 (k < n), -1/2 (k = n), 0 (k > n).

No floating-point cancellation occurs: J has integer entries and D is
half-integer, so the 1D matrices are exact to machine precision.
"""
from __future__ import annotations

import numpy as np
from scipy.special import roots_laguerre, roots_legendre  # noqa: F401  (re-exported)


def jacobi_matrix(size: int) -> np.ndarray:
    """Symmetric Jacobi matrix of orthonormal Laguerre polynomials (multiplication by x)."""
    n = np.arange(size, dtype=float)
    J = np.diag(2.0 * n + 1.0)
    off = -(n[:-1] + 1.0)
    J += np.diag(off, 1) + np.diag(off, -1)
    return J


def derivative_matrix(size: int) -> np.ndarray:
    """D with d/dx [L_n(x) e^{-x/2}] = sum_k D[k, n] L_k(x) e^{-x/2}."""
    D = -np.triu(np.ones((size, size)), 1)
    D -= 0.5 * np.eye(size)
    return D


def moment_matrix(nmax: int, p: int, alpha: float) -> np.ndarray:
    """S_p[a, b] = int_0^inf f_a f_b s^p ds  for a, b <= nmax (exact)."""
    size = nmax + 1 + p
    J = jacobi_matrix(size)
    Jp = np.linalg.matrix_power(J, p) if p > 0 else np.eye(size)
    return Jp[: nmax + 1, : nmax + 1] / alpha ** (p + 1)


def one_d_matrix(nmax: int, p: int, da: int, db: int, alpha: float) -> np.ndarray:
    """int (d^da f_a)(d^db f_b) s^p ds with da, db in {0,1}."""
    S = moment_matrix(nmax, p, alpha)
    D = derivative_matrix(nmax + 1)
    if da:
        S = alpha * (D.T @ S)
    if db:
        S = alpha * (S @ D)
    return S


def evaluate(nmax: int, alpha: float, s: np.ndarray) -> np.ndarray:
    """Values f_n(s) for n = 0..nmax; returns array of shape (nmax+1, len(s)).

    The three-term recurrence is applied directly to f_n = L_n e^{-x/2}
    (it is linear), which is stable for the arguments used here.
    """
    x = alpha * np.asarray(s, dtype=float)
    out = np.empty((nmax + 1,) + x.shape)
    out[0] = np.exp(-0.5 * x)
    if nmax >= 1:
        out[1] = (1.0 - x) * out[0]
    for n in range(1, nmax):
        out[n + 1] = ((2 * n + 1 - x) * out[n] - n * out[n - 1]) / (n + 1)
    return out


def mixed_moment_matrix(nmax_a: int, alpha: float, nmax_b: int, beta: float, p: int) -> np.ndarray:
    """int f_a(s; alpha) f_b(s; beta) s^p ds for different scales (Gauss-Laguerre, exact).

    The integrand is polynomial * exp(-(alpha+beta) s / 2), so a Gauss-Laguerre
    rule with enough nodes integrates it exactly (up to rounding).
    """
    gamma = 0.5 * (alpha + beta)
    nq = (nmax_a + nmax_b + p) // 2 + 2
    x, w = roots_laguerre(nq)
    s = x / gamma
    # f values include exp(-alpha s/2) exp(-beta s/2) = exp(-x); undo it in the weight
    W = w * np.exp(x) / gamma * s ** p
    Fa = evaluate(nmax_a, alpha, s)
    Fb = evaluate(nmax_b, beta, s)
    return (Fa * W) @ Fb.T
