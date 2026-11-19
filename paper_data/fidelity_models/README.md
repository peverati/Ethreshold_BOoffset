# Solvable-model susceptibilities versus M, with the exact Coulomb results
*Last updated: 2026-11-13 20:51*

`fig_models_chi.png`: chi_F (fixed coordinates), chi_inv (global dilation removed), chi_sep (separate dilations of R and rho removed) versus M on 0.05 to 1e4, log-log.
Left panel: `solvers/harmonic_model` Model A (dashed; chi_sep is identically 0 and is not drawn) and Model B with g = 1 (solid). Right panel: exact Coulomb {X+,X+,e-} from `../fidelity/results.csv` (markers), with model B chi_sep (dotted gray) for reference. Horizontal lines: 1/32 (dotted) and 3/64 (dash-dot).

**Convention.** No conversion was needed. `solvers/harmonic_model` (README section 4) and `paper_data/fidelity` both use ln M as variable, F(M, M e^h) = 1 - chi h^2/2, chi = <psi'|psi'> - |<psi|psi'>|^2 with ' = d/d ln M, in fixed Jacobi coordinates. The package columns `chi_fixed`, `chi_inv`, `chi_sep` are used as they are; the Coulomb columns are `chiF`, `chi_inv`, `chi_sep`.

**Data.** `models_chi.csv` (121 log-spaced masses): columns M, A_chiF, A_chi_inv, A_chi_sep, B_g1_chiF, B_g1_chi_inv, B_g1_chi_sep. Generated with the package CLI (run from `solvers/harmonic_model/`): `python3 -m harmonic_model --model A|B --g 1 --masses <grid>`, then `python3 make_fig.py A.csv B.csv` here (the intermediate CSVs were in the scratchpad, not kept; the grid is the first column of `models_chi.csv`). Closed forms, no numerical error beyond double precision.

**Reading the figure.** Large-M limits at M = 1e4: A: chi_F 0.09375 (3/32), chi_inv 0.04687 (3/64), chi_sep 0. B: 0.03136, 0.03074, 0.03073, approaching 1/32. Coulomb: all three approach 1/32 (0.0331, 0.0306, 0.0306), so it follows model B, not A. Differences from model B: Coulomb chi_inv has a shallow maximum near M = 6.9 (B is monotone); Coulomb chi_F is larger at small M (1.22 versus 0.17 for M = 0.05). Small M: chi_inv, chi_sep go to 0 in all cases.
