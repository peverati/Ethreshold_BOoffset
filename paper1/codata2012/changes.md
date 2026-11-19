# CODATA 2012 vs 2018 masses: what changes (2026-11-13)

Masses: old m_p = 1836.15267343, m_mu = 206.768283 (CODATA 2018); new m_p = 1836.15267245, m_mu = 206.7682843 (CODATA 2012, as King, Longford and Cox 2013). The proton shift is -9.8e-10 (5e-13 relative), the muon shift +1.3e-6 (6e-9 relative). Because the shifts are so small, all changes are far below the scientific accuracy of the paper; they matter only for the 11th to 12th digit of energies. Nothing in the paper's conclusions changes.

## Method and convergence
The physical-mass points were recomputed with the unchanged scripts `paper1/mass_scans/scan.py` and `common.py` (copies in this folder with only the two masses changed): exact L=0 (roots 0 and 1, weights (1,1,2) and nested N = 28..40 for both masses, same N_opt = 20), exact L=1 'A' and L=2 'S' (N = 24, 28, 32), BO and BO+ad via `bo_h2plus`. Same settings and the same extrapolation rule as the original scan, so the convergence level is the same. Run times: 1.5 h wall on a loaded machine with two processes. Results: `results.csv` (old, new, difference, error for every quantity; produced by `derive.py`), raw rows `exact_L0.csv`, `exact_rot.csv`, `bo.csv`, logs `log_exact.txt`, `log_rot.txt`, `run_log.txt`, `PROGRESS.md`. `changes_table.md` (from `changes_check.py`) lists every quoted number with its new value at the same printed precision.

Independent check: our new energies agree with King et al. Table II (their CODATA-2012 masses) in all printed digits: Mu2+ -0.5851260982524 (King -0.585126098252403), H2+ -0.5971390630793 (King -0.5971390630793). The first-order mass shift dE/dlnM * dM/M reproduces the recomputed E0 differences (Mu2+ -6.04e-11, H2+ +1.6e-12) exactly, so fixed-basis sequences in the README tables can be shifted by this amount instead of being rerun.

## Headline numbers
Errors are the extrapolation errors of the original method (a heuristic, not a bound); the H2+ L=1 error of 4.5e-11 is the same as in the original scan.

| quantity | old | new | error (new) | where quoted |
|---|---|---|---|---|
| Mu2+ E0 exact | -0.585126098192 | -0.585126098252 | 1.4e-13 | mass_scans README:74; threebody README:267-271 (sequence), 305; paper_draft table (442) ; rotational_table.txt:18 |
| Mu2+ E1 exact (v=1) | -0.558828429024 | -0.558828429151 | 1.6e-13 | threebody README validation (results/validation.txt:70) |
| Mu2+ E_BO | -0.585961513672 | -0.585961513728 | 1e-12 (DVR) | README:74 (printed 10 digits: unchanged) |
| Mu2+ E_BO+ad | -0.585111176455 | -0.585111176516 | 1e-12 | README:74 (unchanged at 10 digits) |
| Mu2+ E-E_BO / E-E_BO+ad | +8.35415480e-4 / -1.49217370e-5 | +8.35415475e-4 / -1.49217362e-5 | 1.4e-13 | unchanged at printed precision |
| H2+ E0 exact | -0.597139063081 | -0.597139063079 | 2.2e-13 | README:76, RESULTS_SUMMARY:48, paper_draft:750/756, threebody README 222-223 |
| H2+ E1 exact (v=1) | -0.5871556791 | -0.587155679096 | 1.8e-12 | |
| H2+ E_BO | -0.597233978551 | -0.597233978550 | 1e-12 | bo_h2plus README:102 |
| H2+ E_BO+ad | -0.597138541555 | -0.597138541554 | 1e-12 | bo_h2plus README:103 |
| H2+ E-E_BO | +9.49154704e-5 | +9.49154705e-5 | 2.2e-13 | |
| H2+ E-E_BO+ad | -5.21525413e-7 | -5.21525315e-7 | 2.2e-13 | |
| Mu2+ vib / d1 / d2 (cm-1) | 5771.67125 / 482.32746 / 1413.96986 | 5771.67123 / 482.32746 / 1413.96985 | 7e-8 / 2e-7 / 7e-8 | not quoted in cm-1 for Mu2+ |
| Mu2+ BO vib / d1 / d2 (cm-1) | 5770.00348 / 482.37502 / 1414.08076 | 5770.00347 / 482.37502 / 1414.08075 | | |
| H2+ vib exact (cm-1) | 2191.09952 | 2191.09952 | 4e-7 | RESULTS_SUMMARY:53 |
| H2+ d1 exact (cm-1) | 58.23195 | 58.23195 (58.231946) | 1e-5 | RESULTS_SUMMARY:53 |
| H2+ d2 exact (cm-1) | 174.2316 | 174.2316 | 1e-6 | |
| H2+ BO vib / d1 / d2 (cm-1) | 2191.14393 / 58.233732 / 174.236889 | 2191.14393 / 58.233732 / 174.236889 | exact (DVR 1e-9 hartree) | bo_h2plus README:106-108 |
| H2+ BO+ad vib (cm-1) | 2191.30003 | 2191.30003 | | bo_h2plus README:106 |
| H2+ <R>, Lindemann (exact) | 2.06391386738, 0.112116591657 | 2.0639138674, 0.112116591673 | | |
| Mu2+ <R>, Lindemann (exact) | 2.2052152336, 0.188946141094 | 2.2052152329, 0.188946140818 | | (diff -7e-10, -3e-10) |
| H2+ d2/d1 | 2.99202751 | 2.99202784 | 5e-7 (L=1 error) | unchanged at 4 decimals |
| Mu2+ d2/d1, M d1 | 2.93155576, 0.4544034 | 2.93155576, 0.4544034 | | |
| M(E-E_BO), H2+ / Mu2+ | 0.174279295 / 0.172737424 | 0.174279295 / 0.172737424 | 4e-10 / 3e-11 | |
| M^1.5(E-E_BO+ad), H2+ / Mu2+ | -0.0410335 / -0.0443655 | -0.0410335 / -0.0443655 | 2e-8 / 4e-10 | |
| BO level counts J=0 | 20 (H2+), 7 (Mu2+) | 20, 7 | | unchanged |

Level-by-level shifts for fixed-basis tables (apply to every row of the sequence, valid to <1e-15 because the basis is converged to 1e-6 or better and dE/dM is smooth): ground state of Mu2+ -6.040e-11 (so threebody README 267-271 become -0.58512601860249, -0.58512609774690, -0.58512609825192, -0.58512609825255), ground state of H2+ +1.619e-12 (so 222-223 and 278-282 become -0.59713906307925 at N=40 etc., see changes_table.md); first excited S root: Mu2+ -1.270e-10, H2+ +3.96e-12. Excited roots 2 and 3 of H2+ (README:284, "-0.57775190441798, -0.56890849870333") were not recomputed; their shifts are of order 1e-12 (expect the same order as E0 and E1: +1.6e-12, +4.0e-12), so only the 12th decimal can move.

## Which quoted numbers change at their printed precision (see changes_table.md for each)
Change:
- Any 12-14 digit energy of Mu2+ (changes at digit 11: -0.585126098192 -> -0.585126098252). Mu2+ E1 likewise (-0.558828429024 -> -0.558828429151).
- 12-14 digit energies of H2+ (digit 12: -0.597139063081 -> -0.597139063079; sequences in threebody README and paper table 422-424 shift by +1.6e-12; the 4 levels list -0.58715567909955 -> -0.58715567909559).
- H2+ E_BO printed to 12 digits (-0.597233978551 -> -0.597233978550) and E_BO+ad (-0.597138541555 -> -0.597138541554): 1e-12 changes, at the level of the DVR error; the 10-digit print in mass_scans README:76 (-0.5972339786) flips at a rounding boundary to -0.5972339785. Safer: quote 9 digits (-0.597233979, -0.597138542), as RESULTS_SUMMARY already does.
- The capped-N table value -0.597139059345 (threebody README:308, paper_draft:445) -> -0.597139059343.
- The rotational table energies (rotational_table.txt:18-19, 12 digits): Mu2+ L=0,1,2 -> -0.585126098252, -0.582928452634, -0.578683577569; H2+ -> -0.597139063045 (N=32 variational), -0.596873738809 (L=1), -0.596345205477 (L=2). Caution: for H2+ L=1 the old table values came from the rotational folder's N=32 variational run and the new ones from the extrapolated scan value (error 4.5e-11), so those differences (up to 4e-11) are convergence-level, not the mass effect (the mass effect itself is ~1e-12). The same holds for README:61 "(0,1) -0.5968737387851 -> -0.5968737388092", "(0,2) -0.5963452054905 -> -0.5963452054773" (differences 2e-11, 1e-11, within the extrapolation error).

Do not change at printed precision (everything else I could find, 70+ entries): all E-E_BO, E-E_BO+ad printed to 4-5 significant figures, all <R>, Lindemann ratios, dR, dE/dlnM, all spacings in hartree (printed 5-6 digits) and in cm-1 (2191.10, 58.232, 2191.144, 2191.300, 58.234, 174.237, 29.156), d2/d1, M d1, vib/d1, E_th (-0.4975935, -0.4997278), M/(M+1), the M(E-E_BO) and M^1.5(E-E_BO+ad) coefficients, 9- to 11-digit energies (H2+ -0.59713906308, -0.597233979, -0.597138542, -0.597139063), and the level counts.

Mass-constant statements that must be reworded (not numbers): paper_draft table of constants (M_p 1836.15267261 "(lit.)", "solvers use 1836.15267343"; muon "solvers use 206.768283"), the H2+ mass-dependence note (paper_draft ~756) and the discussion of -0.597139063079 vs -0.597139063081 (it becomes: with the CODATA-2012 mass the solver gives -0.5971390630793; the published value at 1836.15267247 is -0.597139063079392, the remaining mass difference being 2e-11 relative in dM and ~6e-15 in E), threebody README (CLI examples with 1836.15267343), bo_h2plus README, mass_scans README:44. The sentence in RESULTS_SUMMARY:48 "H2+ -0.597139063081 (reference -0.597139063079 at the slightly lighter m_p; difference explained by the mass)" becomes: our value at the CODATA-2012 mass is -0.5971390630793, agreeing with the reference to 1e-13 (the mass differs from Tang et al. by 2e-8 in m_p, which moves E by 4e-14).

## mass_scans/common.py (step 4)
Not edited. Changing the constants there is not trivial: `scan.py` keys all CSV rows by M (rounded to 9 digits) and `analyze.py` selects the physical masses by value, so every existing row at 206.768283 and 1836.15267343 (energies, exact_rot, bo, descriptors, spacings, gR, density, `data/L0_M*.npz`, figures) would become inconsistent with the constants or be recomputed alongside as duplicates. Recommendation: keep mass_scans as the CODATA-2018 scan (it does not matter for the scaling fits, critical masses or any non-physical grid point), state in its README that the two physical-mass rows are at CODATA 2018, and take all physical-mass table entries for Paper 1 from `paper1/codata2012/results.csv` (or, equivalently, update common.py and rerun `scan.py exact 206.7682843 1836.15267245; scan.py rot ...; scan.py bo ...` into fresh files, which is exactly what this folder did, ~1.5 h). If you want one consistent set, copy the two rows from here into mass_scans/*.csv and replace the two constants; I did not do that to avoid touching shared files.

## Other points
- D2+, T2+ and other deuteron/triton systems appear only in paper1/king2013 (already CODATA 2012); nothing to recompute.
- Rotational validation against Tang et al. (m_p = 1836.15267247) is unchanged.
- No Fortran/ATOM-MOL-nonBO numbers depend on these constants (H2 runs used other masses/notes).
