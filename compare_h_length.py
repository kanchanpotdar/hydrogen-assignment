"""
compare_h_length.py - sensitivity of the PBR heat transfer to the length scale L in Nu = h L / k_r:
    "particle" L = catalyst particle diameter (0.2 mm)  <- what the paper's TEXT specifies
    "tube"     L = tube inner diameter (6 mm)           <- tested as an alternative
The paper text names the particle diameter. Both options were evaluated; the tube diameter reproduces the
published Fig. 4 profiles substantially better and is therefore retained as the Fig. 4-reproduction setting
(default of hyd.py). This is a documented deviation from the paper text, not something the paper states.
Both runs use the corrected Zehner-Schlunder conductivity and the same 853 K / 120 mL/min case.

Each choice is compared with the paper's Fig. 4:
  (a) digitised simulated curves (species RMSE, temperature RMSE),
  (b) experimental outlet points (mean absolute error),
  (c) the temperature dips quoted in the paper text (about 20 K for 57-4Q, 15 K for 25-4Q).

Needs ./data/ (see validation.py) and hyd.py/validation.py in the same folder.
OUTPUT  h_length_comparison.csv, h_length_comparison.png        RUN  python compare_h_length.py
"""
import csv
import dataclasses
import numpy as np
import hyd as h
from validation import IDX, CATS, load_experiment, load_paper_curves, exit_composition

MODES = ["tube", "particle"]
DIP_TEXT = {"57-4Q": 20.0, "25-4Q": 15.0}          # K, paper section 4.2 (HMMC dip is only described as "deep")


def main(out="."):
    import os
    os.makedirs(out, exist_ok=True)
    EXP, DIG = load_experiment(), load_paper_curves()

    # ---- run the model for every (length scale, catalyst) combination --------------------------
    runs = {}
    for m in MODES:
        reactor = dataclasses.replace(h.SMR_REACTOR, h_length=m)       # same geometry, different L
        for c in CATS:
            runs[(m, c)] = (h.run(catalyst=c, Tf=853.0, flow=120.0, reactor=reactor), reactor)

    # ---- metrics ---------------------------------------------------------------------------------
    rows = []
    for (m, c), (res, reactor) in runs.items():
        z, y, T, b = h.profiles(res); n = z <= b["PBR_end"]            # PBR part only
        ex = exit_composition(res)
        sp_rmse = []                                                    # one RMSE per digitised species curve
        for s in ("CH4", "H2", "CO", "CO2", "H2O"):                     # curves vs paper simulation
            if (s, c) in DIG:
                d = DIG[(s, c)]; d = d[d[:, 0] > 0.05]
                sp_rmse.append(np.sqrt(np.mean((np.interp(d[:, 0], z[n], y[IDX[s]][n]) - d[:, 1]) ** 2)))
        t_rmse, n_t = "", 0
        if ("T", c) in DIG:                                             # temperature curve (visible part, z 0.15-5 mm)
            d = DIG[("T", c)]; d = d[(d[:, 0] > 0.15) & (d[:, 0] < 5)]
            n_t = len(d)
            t_rmse = round(float(np.sqrt(np.mean((np.interp(d[:, 0], z, T) - d[:, 1]) ** 2))), 1)
        rows.append(dict(h_length=m, catalyst=c,
                         U_inlet_W_m2K=round(h.U_pbr(res["F0"], 853.0, reactor)[0]),
                         T_dip_K=round(float(T[n].min() - 853.0), 1), T_dip_paper_text_K=DIP_TEXT.get(c, ""),
                         n_species_curves=len(sp_rmse), species_curve_RMSE=round(float(np.mean(sp_rmse)), 5),
                         T_curve_RMSE_K=t_rmse, n_T_points=n_t,
                         outlet_MAE_vs_experiment=round(float(np.mean([abs(ex[s] - EXP[c][s][0]) for s in ex])), 5)))
    with open(os.path.join(out, "h_length_comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"{'L':9s}{'cat':7s}{'U':>6s}{'dT[K]':>8s}{'paper':>7s}{'nSp':>4s}{'spRMSE':>9s}{'T RMSE':>8s}{'nT':>5s}{'outMAE':>9s}")
    for r in rows:
        print(f"{r['h_length']:9s}{r['catalyst']:7s}{r['U_inlet_W_m2K']:6d}{r['T_dip_K']:8.1f}"
              f"{str(r['T_dip_paper_text_K']):>7s}{r['n_species_curves']:4d}{r['species_curve_RMSE']:9.5f}{str(r['T_curve_RMSE_K']):>8s}{r['n_T_points']:5d}"
              f"{r['outlet_MAE_vs_experiment']:9.5f}")

    # ---- plot: tube (solid) vs particle (dashed) vs paper (thick) ---------------------------------
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"HMMC": "tab:red", "57-4Q": "k", "25-4Q": "tab:green"}
    ls = {"tube": "-", "particle": "--"}
    fig, ax = plt.subplots(2, 3, figsize=(14, 7.5)); ax = ax.ravel()
    for a, s in zip(ax, ["CH4", "H2", "CO", "CO2", "H2O", "T"]):
        for c in CATS:
            for m in MODES:
                z, y, T, b = h.profiles(runs[(m, c)][0]); n = z <= b["PBR_end"]
                a.plot(z[n], (T if s == "T" else y[IDX[s]])[n], ls[m], color=col[c], lw=1.4)
            if (s, c) in DIG:
                d = DIG[(s, c)]; a.plot(d[:, 0], d[:, 1], color=col[c], alpha=0.25, lw=6)
        a.set_title(s); a.set_xlabel("PBR z [mm]")
    ax[5].set_xlim(0, 4)                                               # zoom on the inlet temperature dip
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color="k", ls="-", label="L = tube diameter"),
                        Line2D([], [], color="k", ls="--", label="L = particle diameter"),
                        Line2D([], [], color="k", lw=6, alpha=0.25, label="paper Fig. 4 (digitised)")],
               loc="lower center", ncol=3)
    fig.suptitle("PBR heat-transfer length scale vs paper Fig. 4 (853 K, 120 mL/min)")
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(os.path.join(out, "h_length_comparison.png"), dpi=130)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=".")
    main(ap.parse_args().out)
