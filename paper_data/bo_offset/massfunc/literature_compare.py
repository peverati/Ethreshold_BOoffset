"""Compare our m_v(R), m_r(R) (results_fine.csv) with values quoted in
Bunker & Moss, Mol. Phys. 33, 417 (1977); Bunker, McLarnon & Moss, Mol. Phys. 33, 425 (1977);
Kutzelnigg, Mol. Phys. 105, 2627 (2007); Moss, Mol. Phys. 89, 195 (1996);
Schwenke, J. Chem. Phys. 114, 1693 (2001); Jaquet & Kutzelnigg, Chem. Phys. 346, 69 (2008).
Reads results_fine.csv only; no new electronic calculation.  The Moss section solves the 1D BO
vibration-rotation problem on the tabulated E0_exact(R) (finite differences, about 15 s) to predict his Table 4.

Conventions (see LITERATURE_COMPARISON.md):
  ours:        mu_v = (M + m_v)/2, mu_r = (M + m_r)/2   (M = nuclear mass in m_e; leading order in 1/M)
  Kutzelnigg:  1/mu = (2/M)(1 + A/M), 1/mu_rot = (2/M)(1 + B/M)  ->  m_v = -A, m_r = -B
  Bunker-Moss: mu_v = mu/(1+beta), mu_r = mu/(1+alpha)          ->  m_v = -M beta, m_r = -M alpha
Run: python3 paper_data/bo_offset/massfunc/literature_compare.py
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "results_fine.csv"))))
R = np.array([float(r["R"]) for r in rows])
mv = np.array([float(r["m_v"]) for r in rows])
mr = np.array([float(r["m_r"]) for r in rows])
o = np.argsort(R); R, mv, mr = R[o], mv[o], mr[o]

def at(x):
    return float(np.interp(x, R, mv)), float(np.interp(x, R, mr))

U_ME = 1822.888486  # unified mass unit in electron masses (CODATA; value to ~1e-9 relative)

print("== Kutzelnigg 2007, H2+ (A = -m_v, B = -m_r; electron share = m/(1/2)) ==")
kz = [  # (label, R, A or None, B or None)
    ("simple LCAO, R=2",           2.0, -0.20, -0.06),
    ("simple LCAO, R=2.5",         2.5, -0.29, -0.11),
    ("FH (var. alpha), R_e=2",     2.0, -0.31, -0.14),
    ("simple LCAO, A=-0.5 near R=3", 3.0, -0.5, None),
    ("simple LCAO, min A=-0.6 near R=5", 5.0, -0.6, None),
]
print(f"{'case':36s} {'R':>5s} {'-A_K':>7s} {'m_v ours':>9s} {'-B_K':>7s} {'m_r ours':>9s}")
for lab, x, A, B in kz:
    v, r = at(x)
    sA = f"{-A:7.3f}" if A is not None else "      -"
    sB = f"{-B:7.3f}" if B is not None else "      -"
    print(f"{lab:36s} {x:5.2f} {sA} {v:9.4f} {sB} {r:9.4f}")
i = np.argmax(mv)
print(f"ours: m_v max {mv[i]:.5f} at R = {R[i]:.3f}; m_v first reaches 0.5 at R = "
      f"{R[np.argmax(mv >= 0.5)]:.3f}; R_e = 1.99719: m_v, m_r = {at(1.99719332)}")
print(f"ours: participation at R_e: vib {at(1.99719332)[0]/0.5:.1%}, rot {at(1.99719332)[1]/0.5:.1%}")

print("\n== Bunker, McLarnon & Moss 1977, D2/H2 (two electrons; atomic limit m = 1) ==")
dmu_fit = 0.0002538          # mu_v - mu for D2 (u), their eq. (4)
dmu_atom = 1.0070500 - 1.0067757  # mu_atomic - mu (D2), their eqs (5),(6)
print(f"m_v,eff(H2/D2, constant, fitted to v<=20) = 2*dmu*u = {2*dmu_fit*U_ME:.4f} m_e "
      f"(atomic limit check: {2*dmu_atom*U_ME:.4f})")
mu_beta = -0.0002537
print(f"mu*beta = {mu_beta} u  ->  m_v = -2 mu beta = {-2*mu_beta*U_ME:.4f} m_e")
mu_D2 = 1.0067757 * U_ME
alpha = -0.000045
print(f"Unsold estimate alpha(D2) = {alpha}: m_r = -M alpha (as printed) = {-2*mu_D2*alpha:.3f} m_e; "
      f"x4 (factor discussed in LITERATURE_COMPARISON.md) = {-8*mu_D2*alpha:.3f}")
print(f"fitted mu_r - mu = 0.000014 +- 0.000034 u -> m_r = {2*0.000014*U_ME:.3f} +- {2*0.000034*U_ME:.3f} m_e")

print("\n== Bunker & Moss 1977, H2 order-of-magnitude statements ==")
mu_H2 = 0.5036382 * U_ME
print(f"alpha(H2) ~ -1.6e-5 -> m_r = -M alpha = {2*mu_H2*1.6e-5:.4f} (printed formula) or "
      f"{4*2*mu_H2*1.6e-5:.4f} (x4)")
print(f"|beta(H2)| up to ~1e-3 -> m_v up to ~{2*mu_H2*1e-3:.2f}")
print("for comparison H2+ (ours, M = 1836.15): beta(R_e) = -m_v/M = %.2e, alpha(R_e) = -m_r/M = %.2e"
      % (-at(1.99719332)[0]/1836.15267, -at(1.99719332)[1]/1836.15267))

# ---------------------------------------------------------------------------------------------
print("\n== Moss 1996, Mol. Phys. 89, 195: H2+ non-adiabatic corrections (Tables 4, 5, 8) ==")
# Moss gives no mass FUNCTIONS. He fits CONSTANT mu_v, mu_r (and a potential scale 1+gamma) to
# exact non-adiabatic level energies [Moss 1993]. Convention: participation p = (mu - mu_nuc)/(mu_at - mu_nuc),
# mu_at - mu_nuc = 1/2 electron mass / 2  ->  our m = p/2 (atomic limit m = 1/2).
mu_nuc, mu_at = 0.503638235, 0.503775380   # u, Table 5 caption
fits = [  # (label, mu_v [u] or None, mu_r [u] or None, gamma*1e5)
    ("Table 5 model 3 (mu_v only)",          0.503779,  None,      0.0),
    ("Table 5 model 4 (mu_r only)",          None,      0.503719,  0.0),
    ("Table 5 model 6 (mu_v, mu_r)",         0.5037518, 0.5036853, 0.0),
    ("Table 5 model 7 (mu_v, gamma)",        0.503744,  None,      2.38),
    ("Table 5 model 8 (mu_r, gamma)",        None,      0.503671,  4.79),
    ("Table 5 model 9 (all three, 462 lv)",  0.5037686, 0.5036953, -1.51),
    ("Table 6 H2+ N<=5 (all three)",         0.5037604, 0.503676,  -0.627),
]
print(f"{'fit':38s} {'p_v':>6s} {'m_v,eff':>8s} {'p_r':>6s} {'m_r,eff':>8s}")
for lab, v, r, g in fits:
    pv = (v - mu_nuc) / (mu_at - mu_nuc) if v else None
    pr = (r - mu_nuc) / (mu_at - mu_nuc) if r else None
    f = lambda p: (f"{p:6.3f} {p/2:8.4f}" if p is not None else "     -        -")
    print(f"{lab:38s} {f(pv)} {f(pr)}")
print("Table 8 (as printed): p_v = 0.951(8) all N, 0.891(7) N<=5, D2+ 0.884(6); p_r = 0.416(5), 0.28(10), 0.29(15)")
print("Table 7: b = -12.21e-5 u, a = -3.8e-5 u  ->  p_v = -b/(mu_at-mu_nuc) = %.3f, p_r = %.2f"
      % (12.21e-5 / (mu_at - mu_nuc), 3.8e-5 / (mu_at - mu_nuc)))

# --- Our prediction of Moss's Table 4 (first order in m(R)/M; BO vibrational functions) ---------
# H_eff kinetic: (1/M)(1 - m_v(R)/M) |chi'|^2 + N(N+1)(1/M)(1 - m_r(R)/M) chi^2/R^2  (mu = M/2), so
#   Delta E(v,N) = -(1/M^2) [ int m_v chi'^2 dR + N(N+1) int m_r chi^2 / R^2 dR ] .
# Moss's Table 4 = (exact non-adiabatic) - (partitioned adiabatic) dissociation energy = -Delta E.
# Omitted here: the second-order non-adiabatic POTENTIAL Delta W(R) of Bunker-Moss, the O(1/M) R -> mu_e R
# argument scaling, and adiabatic corrections to chi. This is a 1D finite-difference solve of the
# BO potential E0_exact(R) + 1/R from results_fine.csv (h = 0.002, box 0.4-60 bohr; checked vs h = 0.004, box 40); no new electronic calculation.
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
E0 = np.array([float(r["E0_exact"]) for r in rows])[o]
V = CubicSpline(R, E0 + 1.0 / R)
Mv = CubicSpline(R, mv); Mr = CubicSpline(R, mr)
M = 1836.15267343; HCM = 219474.6313705
h = 0.002; x = np.arange(0.4, 60.0 + h / 2, h)
def levels(N):
    d = 2.0 / (M * h * h) + V(x) + N * (N + 1) / (M * x * x)     # -(1/2mu) d2/dR2 = -(1/M) d2/dR2, mu = M/2
    e = np.full(len(x) - 1, -1.0 / (M * h * h))
    w, U = eigh_tridiagonal(d, e, select="v", select_range=(-1.0, -0.5 - 1e-9))
    out = []
    for k in range(len(w)):
        c = U[:, k] / np.sqrt(h)
        dc = np.diff(c) / h; xm = 0.5 * (x[1:] + x[:-1])
        Iv = np.sum(Mv(xm) * dc * dc) * h
        Ir = N * (N + 1) * np.sum(Mr(x) * c * c / x**2) * h
        Tv = np.sum(dc * dc) * h; Tr = N * (N + 1) * np.sum(c * c / x**2) * h
        out.append((k, w[k], Iv * HCM / M**2, Ir * HCM / M**2, Tv * HCM / M**2, Tr * HCM / M**2))
    return out
# Moss Table 4, all 462 printed entries (v: values for N = 0, 1, ... in 1e-3 cm-1), parsed from the PDF text
# and spot-checked against the page.  Two entries look like misprints (v=0,N=2: 0.177; v=10,N=4: 1.037).
MOSS_T4 = {
    0: [109, 111, 177, 123, 133, 145, 158, 174, 191, 211, 232, 258, 283, 311, 342, 373, 406, 441, 478, 516, 555, 596, 638, 682, 727, 772, 819, 867, 916, 966, 1016, 1067, 1118, 1169, 1220, 1271, 1321, 1369, 1414, 1454, 1482],
    1: [311, 311, 314, 322, 330, 342, 354, 367, 383, 402, 422, 443, 468, 493, 520, 548, 578, 610, 643, 677, 713, 750, 789, 828, 869, 910, 952, 995, 1039, 1082, 1126, 1170, 1213, 1255, 1296, 1334, 1368, 1395, 1410],
    2: [492, 494, 496, 504, 510, 521, 530, 544, 559, 575, 593, 612, 634, 656, 680, 706, 733, 762, 792, 822, 854, 888, 921, 957, 992, 1028, 1065, 1102, 1139, 1175, 1210, 1245, 1277, 1306, 1330, 1346, 1349, 1088],
    3: [657, 658, 662, 667, 674, 682, 691, 704, 716, 730, 746, 764, 783, 803, 824, 847, 871, 896, 922, 949, 977, 1006, 1036, 1065, 1095, 1125, 1155, 1184, 1212, 1239, 1263, 1284, 1299, 1306, 1298, 1300],
    4: [806, 807, 810, 814, 820, 827, 836, 846, 857, 870, 883, 898, 915, 932, 951, 971, 991, 1012, 1034, 1057, 1081, 1104, 1127, 1151, 1174, 1197, 1218, 1237, 1254, 1266, 1273, 1271, 1256, 1217],
    5: [938, 938, 940, 944, 949, 956, 963, 971, 981, 991, 1003, 1016, 1029, 1044, 1059, 1075, 1092, 1109, 1126, 1144, 1162, 1179, 1196, 1211, 1226, 1238, 1247, 1252, 1251, 1242, 1219, 1172],
    6: [1051, 1053, 1055, 1057, 1062, 1066, 1072, 1079, 1087, 1095, 1104, 1114, 1125, 1136, 1147, 1159, 1171, 1183, 1195, 1206, 1217, 1227, 1234, 1240, 1242, 1240, 1232, 1215, 1186, 1134],
    7: [1148, 1148, 1149, 1152, 1155, 1159, 1163, 1167, 1173, 1179, 1186, 1192, 1199, 1206, 1213, 1220, 1227, 1232, 1237, 1241, 1242, 1241, 1237, 1228, 1214, 1190, 1154, 1098, 1013],
    8: [1225, 1225, 1226, 1227, 1229, 1231, 1233, 1236, 1238, 1241, 1244, 1247, 1250, 1252, 1254, 1254, 1254, 1252, 1248, 1241, 1230, 1215, 1194, 1165, 1123, 1062, 970],
    9: [1280, 1280, 1280, 1280, 1280, 1281, 1281, 1281, 1280, 1279, 1278, 1277, 1274, 1270, 1265, 1257, 1248, 1236, 1220, 1199, 1171, 1136, 1089, 1024, 926],
    10: [1311, 1311, 1310, 1309, 1037, 1305, 1303, 1299, 1295, 1290, 1284, 1276, 1266, 1255, 1241, 1223, 1202, 1176, 1143, 1102, 1049, 980, 881],
    11: [1315, 1314, 1312, 1309, 1306, 1301, 1295, 1287, 1278, 1267, 1255, 1240, 1221, 1200, 1174, 1143, 1105, 1059, 1002, 928, 826],
    12: [1287, 1285, 1282, 1277, 1271, 1262, 1252, 1240, 1224, 1207, 1186, 1161, 1131, 1096, 1054, 1004, 942, 864, 759],
    13: [1222, 1219, 1214, 1207, 1197, 1185, 1169, 1150, 1127, 1100, 1068, 1030, 985, 931, 866, 784, 675],
    14: [1114, 1110, 1103, 1093, 1079, 1061, 1038, 1011, 978, 938, 891, 835, 767, 682, 568],
    15: [956, 951, 942, 927, 908, 882, 851, 812, 765, 708, 638, 550, 430],
    16: [742, 735, 722, 702, 675, 640, 596, 541, 471, 380],
    17: [466, 458, 440, 413, 376, 326, 260],
    18: [154, 144, 123, 90],
    19: [9, 6],
}
cache = {}
def lev(N):
    if N not in cache:
        cache[N] = levels(N)
    return cache[N]
print("\nTable 4 (N=0): Moss exact non-adiabatic increase vs ours (mass functions only, first order)")
print(f"{'v':>3s} {'E_v (cm-1 rel. limit)':>21s} {'Moss':>7s} {'ours':>7s} {'diff':>7s} {'<m_v>_T':>8s}")
for (k, w, Iv, Ir, Tv, Tr) in lev(0):
    print(f"{k:3d} {(w + 0.5) * HCM:21.2f} {MOSS_T4[k][0]/1e3:7.3f} {Iv:7.3f} {Iv - MOSS_T4[k][0]/1e3:7.3f} {Iv / Tv:8.4f}")
print(f"number of bound N=0 levels in our BO potential: {len(lev(0))} (Moss: v = 0..19)")
diffs, rel, missing, out = [], [], 0, []
for v, vals in MOSS_T4.items():
    for N, val in enumerate(vals):
        L = lev(N)
        if len(L) <= v:            # quasibound in Moss's list; not bound in our 1D box treatment
            missing += 1; continue
        k, w, Iv, Ir, Tv, Tr = L[v]
        d = (Iv + Ir) - val / 1e3
        diffs.append(d); rel.append(d / (val / 1e3))
        if abs(d) > 0.004:
            out.append((v, N, val / 1e3, Iv + Ir))
diffs = np.array(diffs); rel = np.array(rel)
print(f"\nall of Table 4: {len(diffs)} bound levels compared ({missing} quasibound skipped); "
      f"ours - Moss: mean {diffs.mean():+.4f}, rms {np.sqrt(np.mean(diffs**2)):.4f}, "
      f"median |diff| {np.median(abs(diffs)):.4f} cm-1; |diff| <= 0.002 cm-1 for {np.mean(abs(diffs) <= 0.002):.1%}")
bad = np.array([abs(d) > 0.05 for d in diffs])
print(f"excluding the two misprints: max |diff| {abs(diffs[~bad]).max():.4f}, rms {np.sqrt(np.mean(diffs[~bad]**2)):.4f} cm-1, "
      f"max relative {abs(rel[~bad]).max():.1%} (largest at the tiny near-threshold corrections)")
print("entries with |diff| > 0.004 cm-1 (v, N, Moss, ours):", [(v, N, a, round(b, 3)) for v, N, a, b in out])
for vv, Ns in ((0, (1, 3, 5, 10, 20, 30, 35)), (5, (5, 15, 25))):
    print(f"v={vv}: N, Moss, ours (vib + rot), <m_r>_rot")
    for N in Ns:
        k, w, Iv, Ir, Tv, Tr = lev(N)[vv]
        print(f"  {N:3d} {MOSS_T4[vv][N]/1e3:7.3f} {Iv + Ir:7.3f} ({Iv:.3f} + {Ir:.3f}) {Ir / Tr:7.4f}")

# Constant-mass least-squares fit to OUR predicted corrections (analogue of Moss model 6, gamma = 0),
# over all bound levels with N <= 40.  Tells what Moss's procedure would return if Delta W were zero.
A_rows, b_rows = [], []
for N in range(0, 41):
    for (k, w, Iv, Ir, Tv, Tr) in lev(N):
        A_rows.append([Tv, Tr]); b_rows.append(Iv + Ir)
A_ = np.array(A_rows); b_ = np.array(b_rows)
if np.all(A_[:, 1] == 0):
    sol = [b_.sum() / A_[:, 0].sum(), 0.0]
else:
    sol = np.linalg.lstsq(A_, b_, rcond=None)[0]
print(f"\nconstant-mass fit to our corrections ({len(b_)} bound levels, N<=40): m_v,eff = {sol[0]:.4f} "
      f"(p_v = {2*sol[0]:.3f}), m_r,eff = {sol[1]:.4f} (p_r = {2*sol[1]:.3f}); "
      f"rms residual {np.sqrt(np.mean((A_ @ sol - b_)**2)):.3f} cm-1")
A5, b5 = [], []
for N in range(0, 6):
    for (k, w, Iv, Ir, Tv, Tr) in lev(N):
        A5.append([Tv, Tr]); b5.append(Iv + Ir)
A5 = np.array(A5); b5 = np.array(b5)
s5 = np.linalg.lstsq(A5, b5, rcond=None)[0]
# (N<=5 alone cannot separate m_r: the rotational lever arm is tiny, as Moss also found; value not meaningful)
print(f"same, N<=5 ({len(b5)} levels): m_v,eff = {s5[0]:.4f} (p_v = {2*s5[0]:.3f}), m_r,eff = {s5[1]:.4f} (p_r = {2*s5[1]:.3f})")

# ---------------------------------------------------------------------------------------------
print("\n== Schwenke 2001, J. Chem. Phys. 114, 1693: Bunker-Moss functions C_vv, C_rr for H2+ (Fig. 1) ==")
# His effective Hamiltonian, Eq. (15): H_na = d_R^+ C_vv d_R + sum_a N_a^+ C_rr N_a + (C_v cross terms) + C.
# With G_v = -1/(2 mu_R), mu_R = M/2, Eqs. (9),(16),(18):  C_vv = -A_par/mu_R^2 = -m_v/M^2,
# C_rr = -A_perp/(mu_R^2 R^2) = -m_r/(M^2 R^2).  He gives the functions only as a figure (tables in EPAPS,
# not seen). Values below were DIGITIZED from Fig. 1(b) (C_vv x 1e7) and 1(e) (C_rr R^2 x 1e9), 300 dpi
# render, x axis assumed 1..20 bohr between the frame lines; digitization uncertainty about 1% of full scale.
M_moss = 1836.152701   # proton/electron mass ratio used by Moss and (stated) by Schwenke
SCH_CVV = [(1.5, -0.958), (2.0, -1.028), (2.5, -1.131), (3.0, -1.258), (3.5, -1.425), (4.0, -1.589), (5.0, -1.866),
           (6.0, -1.951), (7.0, -1.873), (8.0, -1.751), (9.0, -1.653), (10.0, -1.589), (12.0, -1.528),
           (14.0, -1.507), (16.0, -1.498), (18.0, -1.493), (19.5, -1.493)]
SCH_CRR = [(1.5, -13.96), (2.0, -21.89), (2.5, -31.32), (3.0, -41.76), (3.5, -53.59), (4.0, -66.16), (5.0, -91.32),
           (6.0, -112.70), (7.0, -127.55), (8.0, -136.35), (9.0, -141.38), (10.0, -144.03), (12.0, -146.29),
           (14.0, -147.17), (16.0, -147.55), (18.0, -147.80), (19.5, -147.93)]
print(f"{'R':>5s} {'m_v Schw':>9s} {'m_v ours':>9s} {'m_r Schw':>9s} {'m_r ours':>9s}")
dv, dr = [], []
for (xs, cvv), (_, crr) in zip(SCH_CVV, SCH_CRR):
    sv, sr = -cvv * 1e-7 * M_moss**2, -crr * 1e-9 * M_moss**2
    v, r = at(xs); dv.append(sv - v); dr.append(sr - r)
    print(f"{xs:5.1f} {sv:9.3f} {v:9.4f} {sr:9.4f} {r:9.4f}")
print(f"Schwenke (digitized) - ours: m_v mean {np.mean(dv):+.4f}, max |.| {np.max(np.abs(dv)):.4f}; "
      f"m_r mean {np.mean(dr):+.4f}, max |.| {np.max(np.abs(dr)):.4f}  (asymptote of the digitized curves: "
      f"{-SCH_CVV[-1][1]*1e-7*M_moss**2:.4f}, {-SCH_CRR[-1][1]*1e-9*M_moss**2:.4f}; exact 0.5)")

# ---------------------------------------------------------------------------------------------
print("\n== Jaquet & Kutzelnigg 2008, Chem. Phys. 346, 69: model mass functions, Table 6 ==")
# Their A(R), B(R) are Kutzelnigg's (fixed-trial-function Hylleraas functional) from Finkelstein-Horowitz
# (LCAO, optimized exponent) wave functions; B(pol) adds a p polarization function. m_v = -A, m_r = -B.
JK_T6 = [  # R, -A, -B, -B(pol)  (Jaquet-Kutzelnigg 2008 Table 6; all printed negative)
    (0.2, 0.000657, 0.004421, 0.001997), (0.3, 0.004370, 0.009332, 0.004355), (0.4, 0.011627, 0.015479, 0.007462),
    (0.5, 0.021940, 0.022512, 0.011195), (0.6, 0.034766, 0.030161, 0.015450), (0.7, 0.049585, 0.038199, 0.020120),
    (0.8, 0.066038, 0.046472, 0.025125), (0.9, 0.083832, 0.054871, 0.030402), (1.0, 0.102669, 0.063315, 0.035898),
    (1.1, 0.122309, 0.071758, 0.041575), (1.2, 0.142539, 0.080172, 0.047406), (1.3, 0.163169, 0.088541, 0.053374),
    (1.4, 0.184035, 0.096858, 0.059465), (1.5, 0.204996, 0.105125, 0.065672), (1.6, 0.225931, 0.113344, 0.071992),
    (1.7, 0.246739, 0.121521, 0.078422), (1.8, 0.267332, 0.129663, 0.084963), (1.9, 0.287636, 0.137777, 0.091618),
    (2.0, 0.307589, 0.145868, 0.098387), (2.1, 0.327128, 0.153944, 0.105274), (2.2, 0.346202, 0.162010, 0.112284),
    (2.3, 0.364763, 0.170072, 0.119416), (2.4, 0.382765, 0.178134, 0.126674), (2.5, 0.400182, 0.186196, 0.134062),
    (2.6, 0.416977, 0.194263, 0.141579), (2.7, 0.433111, 0.202335, 0.149225), (2.8, 0.448554, 0.210412, 0.157001),
    (2.9, 0.463278, 0.218491, 0.164905), (3.0, 0.477274, 0.226571, 0.172935), (3.1, 0.490517, 0.234650, 0.181087),
    (3.2, 0.502985, 0.242720, 0.189356), (3.3, 0.514665, 0.250778, 0.197733), (3.4, 0.525546, 0.258816, 0.206211),
    (3.5, 0.535634, 0.266828, 0.214783), (3.6, 0.544930, 0.274806, 0.223437), (3.7, 0.553430, 0.282739, 0.232161),
    (3.8, 0.561143, 0.290619, 0.240941), (3.9, 0.568079, 0.298435, 0.249761), (4.0, 0.574254, 0.306176, 0.258605),
    (4.1, 0.579688, 0.313831, 0.267457), (4.2, 0.584420, 0.321391, 0.276299), (4.3, 0.588466, 0.328843, 0.285112),
    (4.4, 0.591856, 0.336175, 0.293874), (4.5, 0.594621, 0.343376, 0.302569), (4.6, 0.596799, 0.350436, 0.311173),
    (4.7, 0.598427, 0.357344, 0.319667), (4.8, 0.599563, 0.364090, 0.328032), (4.9, 0.600031, 0.370667, 0.336249),
    (5.0, 0.600292, 0.377065, 0.344300), (5.1, 0.600083, 0.383277, 0.352168), (5.2, 0.599500, 0.389297, 0.359836),
    (5.3, 0.598583, 0.395119, 0.367292), (5.4, 0.597372, 0.400740, 0.374521), (5.5, 0.595900, 0.406154, 0.381512),
    (5.6, 0.594199, 0.411361, 0.388257), (5.7, 0.592297, 0.416360, 0.394749), (5.8, 0.590224, 0.421148, 0.400981),
    (5.9, 0.588007, 0.425728, 0.406950), (6.0, 0.585670, 0.430102, 0.412654), (6.1, 0.583235, 0.434272, 0.418092),
    (6.2, 0.580725, 0.438240, 0.423267), (6.3, 0.578158, 0.442012, 0.428181), (6.4, 0.575550, 0.445589, 0.432839),
    (6.5, 0.572919, 0.448981, 0.437244), (6.6, 0.570280, 0.452190, 0.441405), (6.7, 0.567643, 0.455223, 0.445328),
    (6.8, 0.565022, 0.458085, 0.449020), (6.9, 0.562426, 0.460785, 0.452491), (7.0, 0.559864, 0.463326, 0.455749),
    (7.5, 0.547804, 0.473914, 0.469174), (8.0, 0.537369, 0.481555, 0.478661), (8.5, 0.528745, 0.486983, 0.485252),
    (9.0, 0.521856, 0.490794, 0.489780), (9.5, 0.516491, 0.493451, 0.492873), (10.0, 0.512393, 0.495296, 0.494980),
    (11.0, 0.507013, 0.497465, 0.497396), (12.0, 0.504068, 0.498522, 0.498535), (13.0, 0.502479, 0.499058, 0.499092),
    (14.0, 0.501613, 0.499344, 0.499380), (15.0, 0.501128, 0.499510, 0.499542), (16.0, 0.500841, 0.499615, 0.499641),
    (17.0, 0.500659, 0.499687, 0.499709), (18.0, 0.500536, 0.499740, 0.499757), (19.0, 0.500447, 0.499780, 0.499794),
    (20.0, 0.500379, 0.499812, 0.499824), (50.0, 0.500000, 0.500000, 0.500000),
]
T6 = np.array(JK_T6); o6 = np.argsort(T6[:, 0]); T6 = T6[o6]
ov = np.array([at(t)[0] for t in T6[:, 0]]); orr = np.array([at(t)[1] for t in T6[:, 0]])
print(f"{'R':>5s} {'-A FH':>8s} {'m_v ours':>9s} {'-B FH':>8s} {'-B(pol)':>8s} {'m_r ours':>9s}")
for xs in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 15.0, 20.0):
    k = int(np.argmin(abs(T6[:, 0] - xs)))
    print(f"{xs:5.1f} {T6[k,1]:8.4f} {ov[k]:9.4f} {T6[k,2]:8.4f} {T6[k,3]:8.4f} {orr[k]:9.4f}")
kA = int(np.argmax(T6[:, 1]))
print(f"FH vibrational maximum {T6[kA,1]:.4f} at R = {T6[kA,0]} (ours {mv.max():.4f} at {R[np.argmax(mv)]:.3f})")
sel = (T6[:, 0] >= 1.0) & (T6[:, 0] <= 6.0)
print(f"R in [1,6]: max |FH - ours| m_v {np.max(abs(T6[sel,1]-ov[sel])):.3f}, m_r (B) {np.max(abs(T6[sel,2]-orr[sel])):.3f}, "
      f"m_r (B pol) {np.max(abs(T6[sel,3]-orr[sel])):.3f}; ratio B/ours at R=2: {T6[T6[:,0]==2.0,2][0]/at(2.0)[1]:.2f}, "
      f"B(pol)/ours: {T6[T6[:,0]==2.0,3][0]/at(2.0)[1]:.2f}")

# Their level tests: Tables 4 (J=0) and 5 (J=10,20,30) give deviations from Moss's exact corrections (Moss
# 1996 Table 4) for M1 (FH A and B) and M3 (FH A and B(pol)). Signs read from the rendered page (the text
# layer drops minus signs); "positive deviation: correction too large".
JK_M1 = {0: [-.012, -.028, -.033, -.033, -.031, -.026, -.020, -.016, -.013, -.012, -.013, -.017, -.022, -.027,
             -.030, -.030, -.027, -.018, -.006, .010],
         10: [.107, .094, .088, .087, .087, .086, .084, .080, .075, .066, .054, .041, .027, .015, .006, .000],
         20: [.337, .326, .315, .302, .285, .266, .242, .214, .183, .148, .111, .071],
         30: [.515, .483, .445, .397, .340, .271]}
JK_M3 = {10: [.032, .024, .023, .026, .030, .033, .036, .036, .035, .030, .022, .014, .005, -.003, -.008, -.009],
         20: [.118, .123, .127, .129, .128, .124, .116, .103, .088, .070, .050, .031],
         30: [.192, .193, .187, .174, .153, .125]}
# Reproduce their deviations with the same first-order formula, using THEIR functions (spline of Table 6)
# and BO chi from our potential: a check of the convention mapping and of how much of their error is
# due to the model functions.
fA = CubicSpline(T6[:, 0], T6[:, 1]); fB = CubicSpline(T6[:, 0], T6[:, 2]); fBp = CubicSpline(T6[:, 0], T6[:, 3])
def first_order(N, fv, fr):
    out = []
    d = 2.0 / (M * h * h) + V(x) + N * (N + 1) / (M * x * x)
    e = np.full(len(x) - 1, -1.0 / (M * h * h))
    w, U = eigh_tridiagonal(d, e, select="v", select_range=(-1.0, -0.5 - 1e-9))
    xm = 0.5 * (x[1:] + x[:-1])
    clip = lambda f, z: np.where(z < 0.2, 0.0, np.where(z > 50.0, 0.5, f(np.clip(z, 0.2, 50.0))))
    for k in range(len(w)):
        c = U[:, k] / np.sqrt(h); dc = np.diff(c) / h
        out.append((np.sum(clip(fv, xm) * dc * dc) * h + N * (N + 1) * np.sum(clip(fr, x) * c * c / x**2) * h)
                   * HCM / M**2)
    return out
print(f"{'N':>3s} {'v':>3s} {'Moss':>7s} {'FH-Moss':>8s} {'JK M1':>7s} {'FHpol-Moss':>10s} {'JK M3':>7s} {'ours-Moss':>9s}")
agree1, agree3, oursd = [], [], []
for N in (0, 10, 20, 30):
    p1 = first_order(N, fA, fB); p3 = first_order(N, fA, fBp)
    for v, jk in enumerate(JK_M1[N]):
        if v >= len(p1): break
        moss = MOSS_T4[v][N] / 1e3; ours_ = lev(N)[v][2] + lev(N)[v][3]
        d1, d3 = p1[v] - moss, p3[v] - moss
        agree1.append(d1 - jk); oursd.append(ours_ - moss)
        s3 = ""
        if N in JK_M3:
            agree3.append(d3 - JK_M3[N][v]); s3 = f"{d3:10.3f} {JK_M3[N][v]:7.3f}"
        else:
            s3 = f"{d3:10.3f} {'-':>7s}"
        if v % 3 == 0 or N == 30:
            print(f"{N:3d} {v:3d} {moss:7.3f} {d1:8.3f} {jk:7.3f} {s3} {ours_ - moss:9.4f}")
agree1, agree3, oursd = map(np.array, (agree1, agree3, oursd))
print(f"our first-order evaluation of THEIR functions vs their printed deviations: M1 rms {np.sqrt(np.mean(agree1**2)):.4f}, "
      f"max {np.max(abs(agree1)):.4f}; M3 rms {np.sqrt(np.mean(agree3**2)):.4f}, max {np.max(abs(agree3)):.4f} cm-1")
print(f"same {len(oursd)} levels with OUR functions: rms {np.sqrt(np.mean(oursd**2)):.4f}, max {np.max(abs(oursd)):.4f} cm-1")
