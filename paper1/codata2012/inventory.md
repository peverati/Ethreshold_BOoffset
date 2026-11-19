# Inventory: quoted numbers that depend on m_p or m_mu (2026-11-13)

Mass used originally: CODATA 2018 (m_p = 1836.15267343, m_mu = 206.768283) unless stated. New: CODATA 2012 (m_p = 1836.15267245, m_mu = 206.7682843). Line numbers refer to the files as of this date (paper_draft and the student guide were being edited by other agents, so their line numbers may drift; search by the quoted string). Old/new values and the "changes at printed precision?" flag for every entry are in `changes_table.md`.

## Mass constants themselves (text that must be updated)
| file:line | statement | mass |
|---|---|---|
| paper1/mass_scans/common.py:25-26 | `M_PROTON = 1836.15267343`, `M_MUON = 206.768283` | CODATA 2018 |
| paper1/mass_scans/README.md:44 | grid listing "206.768283 ... 1836.15267343" | 2018 |
| solvers/bo_h2plus/README.md:3,77-79 | CLI examples and "M = 1836.15267343 is H2+" | 2018 |
| solvers/threebody/README.md:133,148,222-223,264,273 | CLI examples, Mu2+ (M = 206.768283), H2+ (M = 1836.15267343) | 2018 |
| docs/writeups/paper_draft_emergence.md:423-424, 451 ("M=1836.15267343"), ~727-735 (table: M_p = 1836.15267261 "CODATA (lit.)", "solvers use 1836.15267343", M_mu "solvers use 206.768283"), 756 | | 2018 |
| paper1/RESULTS_SUMMARY.md:25 | "H2+ (M = 1836.15267343)" | 2018 |
| paper1/rotational/log_h2p_L01.txt, validation_h2p_L01/L2.txt, PROGRESS.md:12,28 | m_p = 1836.15267247 (CODATA 2006, Tang et al.): STAYS | 2006 |
| paper1/king2013/* | CODATA 2012 already | 2012 |

## Numbers that depend on m_p / m_mu
All of these come from the physical-mass rows of paper1/mass_scans (exact_L0, exact_rot, bo) or from the solver READMEs' validation runs at the same masses.

### Ground-state energies and BO comparison
| file:line | quantity | old | mass |
|---|---|---|---|
| mass_scans/README.md:74 | Mu2+ E exact, E_BO, E_BO+ad, E-E_BO, E-E_BO+ad | -0.585126098192, -0.5859615137, -0.5851111765, +8.354e-4, -1.492e-5 | 206.768283 |
| mass_scans/README.md:76 | H2+ same five | -0.597139063081, -0.5972339786, -0.5971385416, +9.4915e-5, -5.215e-7 | 1836.15267343 |
| mass_scans/energies.csv, exact_L0.csv, bo.csv, descriptors.csv, spacings.csv, density_cm.csv, gR_*.csv, data/L0_M206.768.npz, L0_M1836.15.npz | the full physical-mass rows | (all columns) | 2018 |
| mass_scans/fits.txt:18-19 | M(E-E_BO) at 206.768, 1836.15 = 0.17274, 0.17428; M^1.5(E-E_BO+ad) = -0.04437, -0.04103 | | 2018 |
| mass_scans/README.md:83-84; paper_draft (577) | 0.17428, -0.0410 (H2+) | | 2018 |
| RESULTS_SUMMARY.md:25 | H2+ BO E0 -0.597233979, BO+ad -0.597138542 | | 2018 |
| RESULTS_SUMMARY.md:28, 48; paper_draft:741, 750, 756 | H2+ -0.59713906308 / -0.597139063081 | | 2018 |
| bo_h2plus/README.md:102-103; paper_draft:457-458 | H2+ BO -0.597233978551; BO+ad -0.597138541555 | | 2018 |
| bo_h2plus/README.md:144 (also mass_scans/README.md:144) | E-E_BO+ad = -5.2e-7 / -5.215e-7 | | 2018 |
| threebody/README.md:222-223, 278-285; results/h2plus.txt; paper_draft:423-428 | H2+ sequence N = 16..40, extrapolations, four levels | -0.59713906308087 etc. | 2018 |
| threebody/README.md:264-271; results/validation.txt:63-73; | Mu2+ sequence N = 12..28, E1 | -0.58512609819215, -0.55882842902454 | 2018 |
| threebody/README.md:305, 308; paper_draft:442, 445 | Mu2+ and H2+ rows of the scan table (E0, E_th) | -0.585126098192, -0.597139059345 (capped-N table, shift applies as below) | 2018 |
| student_guide:256-257, 403 | M/(M+1) = 0.99519, 0.999456; E_th -0.49759, -0.499728 | | 206.77 / 1836.15 (rounded) |

### Spacings, descriptors
| file:line | quantity | old | mass |
|---|---|---|---|
| RESULTS_SUMMARY.md:25 | v0->1 = 2191.14 cm-1, d1 = 58.23, <R> = 2.0639, dR/<R> = 0.1121, 20 levels (Mu2+ 7) | BO | 2018 |
| RESULTS_SUMMARY.md:53; mass_scans/README.md:128 | exact vib 2191.10, d1 58.232 cm-1 | | 2018 |
| mass_scans/README.md:97-98 | <R>, Lindemann (exact/BO/BO+ad), -dE/dlnM at 206.768, 1836.15 | | 2018 |
| mass_scans/README.md:122,124 | vib, E-dependent spacings, M d1, d2/d1, vib/d1 | | 2018 |
| bo_h2plus/README.md:106-110; paper_draft:461-465 | 2191.144, 2191.300 cm-1; 58.234, 174.237; 29.156; 2.06393, 0.23142, 0.1121 | | 2018 |
| rotational/rotational_table.txt:18-19, rotational_table.csv:11-12, scan_L.csv:35-40; threebody/README.md:571-572,578; paper_draft:564-565; student_guide:1415 | L=0,1,2 energies, d1, d2, d2/d1, M d1, dvib, BO d1, BO+ad d1 at 206.768 and 1836.15 | | 2018 |
| student_guide:573, 1382 | "-0.597139", "-0.597234 / -0.597139" (6 digits) | | rounded |

### Level counts (unchanged, listed for completeness)
BO levels J = 0: H2+ 20, Mu2+ 7 (README:25, :130). Depend only on (2M+1)/4 and do not change.

### Not affected
Critical masses (paper1/critical_masses, rotational/critical_L*.txt), mass-scan fits at non-physical masses, all mass_scans rows other than 206.768 and 1836.15, validation against Tang et al. (their mass 1836.15267247), harmonic_model tables (use M = 1836 rounded), paper1/king2013 (already CODATA 2012).
