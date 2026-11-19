"""Assemble the offset table from coupling_results.csv and the exact critical masses."""
import os, pandas as pd, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
c = pd.read_csv(os.path.join(here, "coupling_results.csv"), comment="#")
exact = {(0, 1): 9.525634, (0, 2): 29.63939, (0, 3): 59.31866, (1, 0): 4.2082827, (2, 0): 11.4894217}
MV_INF, MR_E = 0.5, 0.0754            # derived asymptote; m_r(R_e) from rotational spacings (fit, rot_mass_data.csv)
rows = []
for (L, v), Mex in exact.items():
    bo = c[(c.L == L) & (c.c == 0)].reset_index(drop=True)
    ad = c[(c.L == L) & (c.c == 1)].reset_index(drop=True)
    ib = v if L == 0 else 0            # BO list for L=0 includes the v=0 crossing at M=0.965
    ia = v - 1 if L == 0 else 0
    b, a = bo.loc[ib], ad.loc[ia]
    d_ad_1st = -b.Dad / b.B            # first-order BO -> BO+ad shift (derived from D(R))
    lB = a.lam_c * a.B
    frad, fcen = a.Drad / lB, a.Dcen / lB
    off_bo, off_ad = b.M_c - Mex, a.M_c - Mex
    pred_mv = MV_INF * frad
    pred_both_half = MV_INF * (frad + fcen)
    pred_mrRe = MV_INF * frad + MR_E * fcen
    mr_fit = (off_ad - 0.4895 * frad) / fcen if fcen > 0 else np.nan
    rows.append(dict(L=L, v=v, M_exact=Mex, M_BO=b.M_c, M_BOad=a.M_c, off_BO=off_bo, off_BOad=off_ad,
                     shift_ad_exact=b.M_c - a.M_c, shift_ad_firstorder=d_ad_1st, f_rad=frad, f_cen=fcen,
                     pred_nonad_mv_half=pred_mv if L == 0 else pred_mrRe, pred_nonad_all_half=pred_both_half,
                     resid_nonad=off_ad - (pred_mv if L == 0 else pred_mrRe),
                     mv_eff=off_ad / frad if L == 0 else np.nan, mr_eff_given_mv0p4895=mr_fit,
                     pred_off_BO=(b.M_c - a.M_c) + (pred_mv if L == 0 else pred_mrRe)))
df = pd.DataFrame(rows); df.to_csv(os.path.join(here, "offset_table.csv"), index=False)
pd.set_option("display.width", 250); print(df.round(5).to_string())
