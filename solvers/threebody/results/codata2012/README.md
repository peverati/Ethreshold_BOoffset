# threebody results at CODATA 2012
*Last updated: 2026-11-14 09:45*

These are CODATA-2012 copies (m_p = 1836.15267245, m_μ = 206.7682843; the Paper 1 choice, same as King 2013) of `examples/validation_tables.py`, `examples/h2plus_study.py` and `examples/mass_scan.py`. Each writes its output into this folder. The files in `results/` keep the older CODATA-2018 values for the record. The library constant `threebody.system.M_PROTON` is unchanged (still 1836.15267343), so the existing tests are untouched.

Started by `scripts/paper1_day_run.sh`; each script is skipped if its output already exists. Status: done 2026-11-14 08:20–08:29. H₂⁺ ground state −0.59713906307897 (best block, ± 5e−12), which agrees with the Paper 1 value −0.5971390630793 to 3e−13; He −2.90372437703280. The README tables of `solvers/threebody` were already at CODATA 2012 (from `paper1/codata2012/`); these files are the matching saved outputs.
