"""Small-M limit check: as M -> 0 the state tends to a dilated infinite-'nucleus' H^- -like state
(light X pair around the heavy e), psi(M) ~ U_{ln mu} psi_inf, mu = M/(M+1), so
chi_F -> ||D psi_inf||^2 (dln mu/dlnM)^2 = ||D psi_inf||^2/(M+1)^2 and chi_inv, chi_sep -> 0.
Computes ||D psi_inf||^2 for masses (1,1,inf), charges (+1,+1,-1)."""
import numpy as np
import fidcore as fc
from threebody import ThreeBodySystem, optimize_alphas, PerimetricBasis, solve

sysinf = ThreeBodySystem((1.0, 1.0, np.inf), (1.0, 1.0, -1.0), identical=(0, 1), name="Hminus_inf")
al, _, _ = optimize_alphas(sysinf, 14)
out = []
for N in (24, 28):
    b = PerimetricBasis(N, al, "S")
    r = solve(sysinf, b, 1)
    C = b.full_tensor(r.vectors[:, 0])
    q = fc.dilation_quantities(C, C, b)
    out.append(f"H^- (inf mass, N={N}, alphas {al[0]:.5f},{al[2]:.5f}): E = {r.energies[0]:.12f}; ||D psi||^2 = {q['DD']/q['n0']:.10f}")
print("\n".join(out)); open("validation.txt", "a").write("\n".join(out) + "\n")
