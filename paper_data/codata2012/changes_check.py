"""Builds changes_table.md: every quoted old value -> new value (CODATA 2012) and whether it changes at printed precision."""
import csv, os, re, math
HERE = os.path.dirname(os.path.abspath(__file__))
R = {(r["mass"], r["quantity"]): r for r in csv.DictReader(open(os.path.join(HERE, "results.csv")))}
def new(m, q): return float(R[(m, q)]["new"])
def err(m, q): return float(R[(m, q)]["err_new"])
sh_mu, sh_p = new("mu","E0_exact")-float(R[("mu","E0_exact")]["old"]), new("p","E0_exact")-float(R[("p","E0_exact")]["old"])
# E_th = -M/(2(M+1))
eth = lambda M: -M/(2*(M+1))
Eth_mu, Eth_p = eth(206.7682843), eth(1836.15267245)
pat = re.compile(r'^([+\-]?)(\d*\.?\d+)(?:e([+\-]?\d+))?$')
def fmt_like(s, x):
    s0 = s.replace("−", "-").replace("\\", "").replace("*", "")
    m = pat.match(s0)
    if not m: return None
    dec = len(m.group(2).split(".")[1]) if "." in m.group(2) else 0
    if m.group(3) is not None:
        e = int(m.group(3)); v = x / 10 ** e
        out = f"{v:.{dec}f}"
        sign = "+" if m.group(1) == "+" and v >= 0 else ""
        return sign + out + f"e{e}"
    out = f"{x:.{dec}f}"
    return ("+" if m.group(1) == "+" and x >= 0 else "") + out
def norm(s):
    s = s.replace("−", "-").replace("\\", "").replace("*", "")
    m = pat.match(s)
    if m and m.group(3) is not None:
        return f"{m.group(1)}{m.group(2)}e{int(m.group(3))}"
    return s
E = []
def add(loc, q, old, newval, unit=""): E.append((loc, q, old, newval, unit))
def sh(old, shift): return float(old.replace("−", "-")) + shift
P, MU = "p", "mu"
# --- exact energies, mass_scans README / fits
add("mass_scans/README.md:74", "Mu2+ E exact", "−0.585126098192", new(MU,"E0_exact"))
add("mass_scans/README.md:74", "Mu2+ E_BO", "−0.5859615137", new(MU,"E0_BO"))
add("mass_scans/README.md:74", "Mu2+ E_BO+ad", "−0.5851111765", new(MU,"E0_BOad"))
add("mass_scans/README.md:74", "Mu2+ E-E_BO", "+8.354e−4", new(MU,"E0_exact_minus_BO"))
add("mass_scans/README.md:74", "Mu2+ E-E_BOad", "−1.492e−5", new(MU,"E0_exact_minus_BOad"))
add("mass_scans/README.md:76; RESULTS_SUMMARY:48; paper_draft:750,756", "H2+ E exact (1836.15267343)", "−0.597139063081", new(P,"E0_exact"))
add("RESULTS_SUMMARY:28; paper_draft:741", "H2+ E exact (11 dec.)", "−0.59713906308", new(P,"E0_exact"))
add("mass_scans/README.md:76", "H2+ E_BO", "−0.5972339786", new(P,"E0_BO"))
add("mass_scans/README.md:76", "H2+ E_BO+ad", "−0.5971385416", new(P,"E0_BOad"))
add("mass_scans/README.md:76", "H2+ E-E_BO", "+9.4915e−5", new(P,"E0_exact_minus_BO"))
add("mass_scans/README.md:76; bo_h2plus README:144", "H2+ E-E_BOad", "−5.215e−7", new(P,"E0_exact_minus_BOad"))
add("RESULTS_SUMMARY:25; bo_h2plus README:102; paper_draft:457", "H2+ E_BO (12/9 dec.)", "−0.597233978551", new(P,"E0_BO"))
add("RESULTS_SUMMARY:25", "H2+ E_BO (9 dec.)", "−0.597233979", new(P,"E0_BO"))
add("bo_h2plus README:103; paper_draft:458", "H2+ E_BO+ad", "−0.597138541555", new(P,"E0_BOad"))
add("RESULTS_SUMMARY:25", "H2+ E_BO+ad (9 dec.)", "−0.597138542", new(P,"E0_BOad"))
add("fits.txt:18; mass_scans/README.md:83; RESULTS_SUMMARY", "M(E-E_BO), Mu2+", "0.17274", new(MU,"M*(E-E_BO)"))
add("fits.txt:18; mass_scans/README.md:83", "M(E-E_BO), H2+", "0.17428", new(P,"M*(E-E_BO)"))
add("fits.txt:19", "M^1.5(E-E_BOad), Mu2+", "-0.04437", new(MU,"M^1.5*(E-E_BOad)"))
add("fits.txt:19; README:84; paper_draft:577", "M^1.5(E-E_BOad), H2+", "-0.04103", new(P,"M^1.5*(E-E_BOad)"))
add("README:84; paper_draft:577", "M^1.5(E-E_BOad), H2+ (3 dec.)", "-0.0410", new(P,"M^1.5*(E-E_BOad)"))
# descriptors
add("mass_scans/README.md:97", "Mu2+ <R>", "2.20522", new(MU,"exact_R"))
add("mass_scans/README.md:97", "Mu2+ Lindemann exact", "0.18895", new(MU,"exact_lindemann"))
add("mass_scans/README.md:97", "Mu2+ Lindemann BO", "0.18908", new(MU,"BO_lindemann"))
add("mass_scans/README.md:97", "Mu2+ Lindemann BO+ad", "0.18902", new(MU,"BOad_lindemann"))
add("mass_scans/README.md:98", "H2+ <R> exact", "2.06391", new(P,"exact_R"))
add("mass_scans/README.md:98", "H2+ Lindemann exact", "0.11212", new(P,"exact_lindemann"))
add("mass_scans/README.md:98", "H2+ Lindemann BO", "0.11212", new(P,"BO_lindemann"))
add("RESULTS_SUMMARY:25; bo_h2plus README:110", "H2+ <R> (BO, 5 dec.)", "2.06393", new(P,"BO_R"))
add("RESULTS_SUMMARY:25", "H2+ <R> (4 dec.)", "2.0639", new(P,"BO_R"))
add("bo_h2plus README:110; paper_draft:465", "H2+ dR (BO)", "0.23142", new(P,"BO_dR"))
add("RESULTS_SUMMARY:25; bo_h2plus README:110", "H2+ dR/<R> (BO)", "0.1121", new(P,"BO_lindemann"))
add("mass_scans/README.md:98", "H2+ -dE/dlnM", "0.002847", -new(P,"exact_dE_dlnM"))
add("mass_scans/README.md:97", "Mu2+ -dE/dlnM", "0.009579", -new(MU,"exact_dE_dlnM"))
# spacings
add("mass_scans/README.md:122", "Mu2+ vib (hartree)", "2.62977e−2", new(MU,"vib_hartree"))
add("mass_scans/README.md:122; paper_draft", "Mu2+ d1", "2.19765e−3", new(MU,"d1_hartree"))
add("mass_scans/README.md:122", "Mu2+ M d1", "0.4544", new(MU,"M*d1"))
add("mass_scans/README.md:122; rotational_table.txt:18; paper_draft:564", "Mu2+ d2/d1", "2.9316", new(MU,"d2/d1"))
add("mass_scans/README.md:122", "Mu2+ vib/d1", "11.97", new(MU,"vib/d1"))
add("mass_scans/README.md:122", "Mu2+ rigid 2<1/R2> (as M d1 col.)", "0.460", new(MU,"rigid_rotor_2<1/R2>/M_hartree")*206.7682843)
add("mass_scans/README.md:124", "H2+ vib (hartree)", "9.98338e−3", new(P,"vib_hartree"))
add("mass_scans/README.md:124; rotational_table.txt:19", "H2+ d1 (hartree)", "2.65324e−4", new(P,"d1_hartree"))
add("mass_scans/README.md:124", "H2+ M d1", "0.4872", new(P,"M*d1"))
add("mass_scans/README.md:124; threebody README:578", "H2+ rigid 2<1/R2> M", "0.4878", new(P,"rigid_rotor_2<1/R2>/M_hartree")*1836.15267245)
add("mass_scans/README.md:124; paper_draft:565; student:1415", "H2+ d2/d1", "2.9920", new(P,"d2/d1"))
add("mass_scans/README.md:124", "H2+ vib/d1", "37.6", new(P,"vib/d1"))
add("mass_scans/README.md:128; RESULTS_SUMMARY:53", "H2+ vib exact (cm-1)", "2191.10", new(P,"vib_cm-1"))
add("mass_scans/README.md:128; RESULTS_SUMMARY:53", "H2+ d1 exact (cm-1)", "58.232", new(P,"d1_cm-1"))
add("RESULTS_SUMMARY:25", "H2+ v 0->1 (BO, 2 dec.)", "2191.14", new(P,"BO_vib_cm-1"))
add("RESULTS_SUMMARY:25", "H2+ d1 (BO, 2 dec.)", "58.23", new(P,"BO_d1_cm-1"))
add("bo_h2plus README:106; paper_draft:461", "H2+ v 0->1 BO", "2191.144", new(P,"BO_vib_cm-1"))
add("bo_h2plus README:106; paper_draft:461", "H2+ v 0->1 BO+ad", "2191.300", new(P,"BOad_vib_cm-1"))
add("bo_h2plus README:108; paper_draft:463", "H2+ d1 BO", "58.234", new(P,"BO_d1_cm-1"))
add("bo_h2plus README:108; paper_draft:463", "H2+ d2 BO", "174.237", new(P,"BO_d2_cm-1"))
add("bo_h2plus README:109; paper_draft:464", "H2+ <1/R2>/(2muN) BO", "29.156", new(P,"BO_<1/R2>/(2muN)_cm-1"))
add("rotational_table.txt:18 (d2)", "Mu2+ d2 (hartree)", "6.44252e-03", new(MU,"d2_hartree"))
add("rotational_table.txt:19 (d2)", "H2+ d2 (hartree)", "7.93858e-04", new(P,"d2_hartree"))
add("rotational_table.txt:18", "Mu2+ M d1", "0.45440", new(MU,"M*d1"))
add("rotational_table.txt:19", "H2+ M d1", "0.48718", new(P,"M*d1"))
add("rotational_table.txt:19", "H2+ d2/d1", "2.9920", new(P,"d2/d1"))
add("rotational_table.txt:18 (dvib)", "Mu2+ dvib", "2.6298e-02", new(MU,"vib_hartree"))
add("rotational_table.txt:19 (dvib)", "H2+ dvib", "9.9834e-03", new(P,"vib_hartree"))
add("rotational_table.txt:18 (BO d1)", "Mu2+ BO d1", "2.19786e-03", new(MU,"BO_d1_hartree"))
add("rotational_table.txt:19 (BO d1)", "H2+ BO d1", "2.65332e-04", new(P,"BO_d1_hartree"))
add("rotational_table.txt:18 (BO+ad d1)", "Mu2+ BOad d1", "2.19839e-03", new(MU,"BOad_d1_hartree"))
add("rotational_table.txt:19 (BO+ad d1)", "H2+ BOad d1", "2.65335e-04", new(P,"BOad_d1_hartree"))
add("threebody README:571; paper_draft", "Mu2+ M d1 (3 dec.)", "0.4544", new(MU,"M*d1"))
# thresholds
add("threebody README:305,571; paper_draft:442", "Mu2+ E_th", "−0.4975935", Eth_mu)
add("threebody README:308,572; paper_draft:445; student:257", "H2+ E_th", "−0.4997278", Eth_p)
add("student:256", "Mu2+ M/(M+1)", "0.99519", 206.7682843/207.7682843)
add("student:257", "H2+ M/(M+1)", "0.999456", 1836.15267245/1837.15267245)
# threebody README sequences at fixed basis (first-order mass shift, validated by actual E0 shifts)
for old in ["−0.58512601854209", "−0.58512609768650", "−0.58512609819152", "−0.58512609819215"]:
    add("threebody README:267-271; validation.txt:63-67", "Mu2+ sequence (fixed basis), shift %.3e" % sh_mu, old, sh(old, sh_mu))
add("threebody README:305; paper_draft:442", "Mu2+ extrapolated", "−0.585126098192", new(MU,"E0_exact"))
for old in ["−0.59713244875271","−0.59713869065408","−0.59713904797634","−0.59713906246867","−0.59713906306281","−0.59713905530895","−0.59713906291303","−0.59713906307571","−0.59713906308049","−0.59713906308087"]:
    add("threebody README:222-223,278-282; paper_draft:423-424; results/h2plus.txt", "H2+ sequence (fixed basis), shift %.3e" % sh_p, old, sh(old, sh_p))
add("threebody README:222; paper_draft:423", "H2+ extrapolated (3417)", "−0.5971390630882", new(P,"E0_exact") + 0.5971390630882-0.59713906308091 if False else sh("−0.5971390630882", sh_p))
add("threebody README:223; paper_draft:424", "H2+ extrapolated (3311)", "−0.59713906308091", new(P,"E0_exact"))
add("threebody README:284; paper_draft:428", "H2+ E1 (N=40)", "−0.58715567909955", sh("−0.58715567909955", new(P,"E1_exact")-float(R[(P,"E1_exact")]["old"])))
add("threebody README:243", "H2+ v=1", "−0.587155679", new(P,"E1_exact"))
add("threebody README:305", "Mu2+ E1 (not quoted at 12 digits in README; validation.txt:70)", "−0.55882842902454", new(MU,"E1_exact"))
add("threebody README:308; paper_draft:445", "H2+ capped-table E0 (err 4e-8)", "−0.597139059345", sh("−0.597139059345", sh_p))
add("rotational_table.txt:18", "Mu2+ L=0 (N=32 variational)", "−0.585126098192", sh("−0.585126098192", sh_mu))
add("rotational_table.txt:19", "H2+ L=0 (N=32 variational)", "−0.597139063047", sh("−0.597139063047", sh_p))
add("rotational_table.txt:18", "Mu2+ L=1", "−0.582928452560", new(MU,"E0_exact")*0 + float(open(os.path.join(HERE,"exact_rot.csv")).read().split("\n")[1].split(",")[0]) * 0 + 0 if False else -0.5829284526336)
add("rotational_table.txt:18", "Mu2+ L=2", "−0.578683577472", -0.5786835775695)
add("rotational_table.txt:19", "H2+ L=1", "−0.596873738766", -0.5968737388092)
add("rotational_table.txt:19", "H2+ L=2", "−0.596345205474", -0.5963452054773)
add("mass_scans/README.md:61 (L=1; physical-mass scan value)", "H2+ (0,1) at 1836.15267343", "−0.5968737387851", -0.5968737388092)
add("mass_scans/README.md:61", "H2+ (0,2) at 1836.15267343", "−0.5963452054905", -0.5963452054773)
rows = []
for loc, q, old, nv, unit in E:
    f = fmt_like(old, nv)
    ch = (f is not None) and norm(f) != norm(old)
    rows.append((loc, q, old, f if f else repr(nv), "YES" if ch else "no"))
with open(os.path.join(HERE, "changes_table.md"), "w") as fh:
    fh.write("| where quoted | quantity | old printed | new at same precision | changes? |\n|---|---|---|---|---|\n")
    for r in rows: fh.write("| " + " | ".join(r) + " |\n")
print(open(os.path.join(HERE, "changes_table.md")).read())
