"""Symbolic (sympy) checks of the derivations."""
import sympy as sp

M, k, kp, g, R = sp.symbols("M k kp g R", positive=True)


def test_jacobi_kinetic_and_springs():
    # rows: R_cm, R, rho in terms of r1, r2, r3 (one Cartesian component)
    Mt = 2 * M + 1
    J = sp.Matrix([[M / Mt, M / Mt, 1 / Mt], [-1, 1, 0], [-sp.Rational(1, 2), -sp.Rational(1, 2), 1]])
    Minv = sp.diag(1 / M, 1 / M, 1)
    T = sp.simplify(J * Minv * J.T)          # T = (1/2) p^T T p in new coordinates
    assert sp.simplify(T - sp.diag(1 / Mt, 2 / M, (2 * M + 1) / (2 * M))) == sp.zeros(3)
    # springs
    Rv, rv = sp.symbols("Rv rv")
    r1 = -Rv / 2 - rv / Mt * 0  # midpoint frame: r1 = m - R/2, r2 = m + R/2, r3 = m + rho
    r1, r2, r3 = -Rv / 2, Rv / 2, rv
    assert sp.expand((r3 - r1) ** 2 + (r3 - r2) ** 2 - (2 * rv ** 2 + Rv ** 2 / 2)) == 0


def test_cm_frame_coordinates():
    x1, x2, x3, Rv, rv = sp.symbols("x1 x2 x3 Rv rv")
    # x_i = r_i - R_cm with M x1 + M x2 + x3 = 0
    sol = sp.solve([sp.Eq(Rv, x2 - x1), sp.Eq(rv, x3 - (x1 + x2) / 2), sp.Eq(M * x1 + M * x2 + x3, 0)],
                   [x1, x2, x3], dict=True)[0]
    assert sp.simplify(sol[x1] - (-Rv / 2 - rv / (2 * M + 1))) == 0
    assert sp.simplify(sol[x3] - 2 * M * rv / (2 * M + 1)) == 0
    # x1|x3 parametrization
    s = sp.solve([sp.Eq(M * x1 + M * x2 + x3, 0)], [x2], dict=True)[0]
    Rexpr = sp.simplify((x2 - x1).subs(s))
    rexpr = sp.simplify((x3 - (x1 + x2) / 2).subs(s))
    assert sp.simplify(Rexpr - (-2 * x1 - x3 / M)) == 0
    assert sp.simplify(rexpr - (2 * M + 1) * x3 / (2 * M)) == 0


def test_model_b_radial_solution():
    K = k + 2 * kp
    lam = (sp.sqrt(1 + 4 * g * M) - 1) / 2
    a = sp.sqrt(M * K) / 2
    u = R ** (lam + 1) * sp.exp(-a * R ** 2 / 2)
    E = sp.sqrt(K / M) * (lam + sp.Rational(3, 2))
    res = -sp.diff(u, R, 2) / M + (g / R ** 2 + K * R ** 2 / 4) * u - E * u
    assert sp.simplify(sp.expand(res / u)) == 0


def test_energy_derivatives_symbolic():
    from harmonic_model import model_b as B
    t = sp.symbols("t", real=True)
    Mt = sp.exp(t)
    K = 1 + 2 * sp.Rational(3, 10)
    lam = (sp.sqrt(1 + 4 * sp.Rational(7, 10) * Mt) - 1) / 2
    E = sp.sqrt(K / Mt) * (lam + sp.Rational(3, 2)) + sp.Rational(3, 2) * sp.sqrt((2 * Mt + 1) / Mt)
    for MM in (0.1, 1.0, 37.0):
        tv = sp.log(sp.Float(MM, 30))
        d1 = float(sp.diff(E, t).subs(t, tv))
        d2 = float(sp.diff(E, t, 2).subs(t, tv))
        assert abs(d1 - B.dE_dlnM(MM, 0.7, 1.0, 0.3)) < 1e-12
        assert abs(d2 - B.d2E_dlnM2(MM, 0.7, 1.0, 0.3)) < 1e-12
