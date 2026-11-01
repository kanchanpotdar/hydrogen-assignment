"""
hyd_contours.py - operating-window and axial contour figures for the SMR model (hyd.py).

OUTPUT FIGURES
  contours_operating_window.png   X_CH4, outlet y_H2 and CO/converted-CH4 ratio versus furnace
                                  temperature (700-1100 K) and feed flow (60-200 mL/min), 3 catalysts.
  contours_axial_<catalyst>.png   y_CH4, y_H2, y_CO, y_CO2, y_H2O and T - T_furnace on a map of
                                  PBR position z (0-12 mm) versus furnace temperature, at 120 mL/min.

IMPORTANT NOTES (state them in the report)
  * The model is 1-D (axial only). The axial "contours" are z-versus-furnace-temperature maps,
    NOT two-dimensional spatial (r, z) contours.
  * Operating-window values are taken at the physical P-R2 outlet (40 mm); the paper's Fig. 4 shows only
    the first 3 mm of P-R2.
  * Supplied Figure 4 comparisons cover 853 K and 120 mL/min only. Red stars mark
    paper experimental conditions, including 80 mL/min without supplied validation data here.
  * "CO / converted CH4" = F_CO,out / (F_CH4,in - F_CH4,out). It is this work's own metric, NOT a quantity
    of the paper and not called "selectivity". It is masked where X_CH4 < 5 % (ratio undefined).
  * A failed simulation point is reported with its conditions and the error is re-raised
    (no silent NaN contours).

RUN   python hyd_contours.py [--out folder] [--workers N]   (needs hyd.py in the same folder)
      Runtime depends on the hardware and N: about 3 minutes with 2 workers on a 2-core machine; allow several
      minutes. Default N = 4 workers (never more than the number of CPU cores).
"""
import os
import sys
import numpy as np
from multiprocessing import Pool
import hyd

# ---- grids and settings ----------------------------------------------------------------------------
CATS = ["HMMC", "57-4Q", "25-4Q"]
T_GRID = np.linspace(700, 1100, 17)         # furnace temperatures of the operating window [K]
FLOW_GRID = np.linspace(60, 200, 11)        # feed flows of the operating window [mL/min at RTP]
FLOW_AXIAL = 120.0                          # flow of the axial maps (the paper's main case)
T_AXIAL = np.linspace(700, 1100, 33)        # furnace temperatures of the axial maps [K]
Z_MM = np.linspace(0.0, 12.0, 241)          # PBR positions of the axial maps [mm] (PBR is 12 mm long)
X_MIN_CO = 0.05                             # mask the CO ratio below this CH4 conversion
PAPER_POINTS = [(853.0, 80.0), (853.0, 120.0)]   # the paper's two experimental conditions


def _point(args):
    """one operating-window point: returns (X_CH4, y_H2 at outlet, CO/converted CH4)"""
    cat, Tf, flow = args
    try:
        res = hyd.run(catalyst=cat, Tf=Tf, flow=flow, n=200)    # full PBR -> P-R1 -> P-R2 chain
        F0 = res["F0"]                                          # inlet molar flows [mol/s]
        Fe = res["PR2"][1][:6, -1]                              # outlet molar flows of the 40 mm P-R2
        conv = (F0[0] - Fe[0]) / F0[0]                          # methane conversion
        yH2 = Fe[4] / Fe.sum()                                  # outlet H2 mole fraction
        ratio = Fe[2] / (F0[0] - Fe[0]) if conv >= X_MIN_CO else np.nan   # CO per CH4 converted
        return conv, yH2, ratio
    except Exception as e:                                      # report WHY and where, then stop
        print(f"FAILED operating point: catalyst={cat} T={Tf:.0f} K flow={flow:.0f} mL/min: {e}", file=sys.stderr)
        raise


def _axial(args):
    """one furnace temperature of the axial maps: returns the PBR profiles resampled on Z_MM"""
    cat, Tf = args
    try:
        res = hyd.run(catalyst=cat, Tf=Tf, flow=FLOW_AXIAL, n=400)
        z, Y = res["PBR"]; z = z * 1e3                          # PBR positions [mm]
        y = hyd.mole_fractions(Y)                               # rows: CH4, H2O, CO, CO2, H2, He
        f = lambda v: np.interp(Z_MM, z, v)                     # resample onto the common grid
        return [f(y[0]), f(y[4]), f(y[2]), f(y[3]), f(y[1]), f(Y[6]) - Tf]   # CH4,H2,CO,CO2,H2O,dT
    except Exception as e:
        print(f"FAILED axial point: catalyst={cat} T={Tf:.0f} K: {e}", file=sys.stderr)
        raise


def fig_operating(pool, out):
    """3x3 grid: rows = catalysts, columns = X_CH4, y_H2, CO ratio; filled + labelled contour lines"""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    jobs = [(c, T, q) for c in CATS for T in T_GRID for q in FLOW_GRID]
    r = np.array(pool.map(_point, jobs, chunksize=4)).reshape(len(CATS), len(T_GRID), len(FLOW_GRID), 3)
    titles = [("Methane conversion X$_{CH4}$ [-]", 0, (0, 1)),            # clipped to [0,1]: it cannot exceed 1
              ("H$_2$ mole fraction at outlet [-]", 1, None),
              ("CO / converted CH$_4$ [-]  (own metric)", 2, None)]
    fig, ax = plt.subplots(3, 3, figsize=(15, 12), constrained_layout=True)
    TT, QQ = np.meshgrid(T_GRID, FLOW_GRID, indexing="ij")
    for i, cat in enumerate(CATS):
        for j, (ttl, k, lim) in enumerate(titles):
            a = ax[i, j]
            Z = np.ma.masked_invalid(r[i, :, :, k])
            if lim is not None:
                Z = np.clip(Z, *lim)                                      # removes tiny numerical overshoot above 1
            lv = np.linspace(Z.min(), Z.max(), 16) if lim is None else np.linspace(lim[0], lim[1], 21)
            cf = a.contourf(TT, QQ, Z, levels=lv, cmap="viridis")         # filled contours
            cs = a.contour(TT, QQ, Z, levels=lv[::2], colors="k", linewidths=0.5)   # labelled lines
            a.clabel(cs, fontsize=6, fmt="%.3g")
            a.plot(*zip(*PAPER_POINTS), "r*", ms=13, mec="w", label="paper experiments")   # paper conditions, not additional validation
            fig.colorbar(cf, ax=a, format="%.3g")
            a.set_title(f"{cat}: {ttl}", fontsize=10)
            a.set_xlabel("furnace temperature [K]"); a.set_ylabel("feed flow [mL/min RTP]")
    ax[0, 0].legend(loc="lower right", fontsize=8)
    fig.suptitle("SMR predictions: full 40 mm P-R2 outlet; stars mark paper conditions.\nSupplied benchmark: 853 K, 120 mL/min only.")
    fig.savefig(os.path.join(out, "contours_operating_window.png"), dpi=130, bbox_inches="tight", pad_inches=.2); plt.close(fig)


def fig_axial(pool, out):
    """one figure per catalyst: 2x3 maps of z versus furnace temperature (CH4, H2, CO, CO2, H2O, T - Tf)"""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    names = [("y$_{CH4}$ [-]", "viridis"), ("y$_{H2}$ [-]", "viridis"), ("y$_{CO}$ [-]", "viridis"),
             ("y$_{CO2}$ [-]", "viridis"), ("y$_{H2O}$ [-]", "viridis"), ("T - T$_{furnace}$ [K]", "RdBu_r")]
    ZZ, TT = np.meshgrid(Z_MM, T_AXIAL)
    for cat in CATS:
        r = pool.map(_axial, [(cat, T) for T in T_AXIAL])                 # one profile set per furnace temperature
        fields = [np.ma.masked_invalid(np.array([x[k] for x in r])) for k in range(6)]
        fig, ax = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True); ax = ax.ravel()
        for a, F, (ttl, cmap) in zip(ax, fields, names):
            if cmap == "RdBu_r":                                          # temperature deviation: asymmetric scale
                lo, hi = min(F.min(), -1.0), max(F.max(), 1.0)            # so the cold dip is not hidden
                norm = TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=hi)
                cf = a.contourf(ZZ, TT, F, levels=np.linspace(lo, hi, 25), cmap=cmap, norm=norm)
                a.contour(ZZ, TT, F, levels=[-1.0], colors="k", linewidths=0.8)   # -1 K line = thermally recovered
            else:
                cf = a.contourf(ZZ, TT, F, levels=20, cmap=cmap)
                a.contour(ZZ, TT, F, levels=8, colors="k", linewidths=0.4)
            fig.colorbar(cf, ax=a, format="%.3g"); a.set_title(f"{cat}: {ttl}")
            a.set_xlabel("PBR axial position z [mm]"); a.set_ylabel("furnace temperature [K]")
            a.axhline(853, c="r", ls="--", lw=0.8)                         # 853 K = the validated case
        fig.suptitle(f"{cat}, {FLOW_AXIAL:.0f} mL/min: axial maps vs furnace temperature (red dashed = 853 K; black contour = -1 K recovery). 1-D model.")
        fig.savefig(os.path.join(out, f"contours_axial_{cat}.png"), dpi=120, bbox_inches="tight", pad_inches=.2); plt.close(fig)


def main(out=None, workers=4):
    import multiprocessing as mp
    from run_metadata import prepare_output, record_run
    from dataclasses import asdict
    out = prepare_output(out, "contours")
    workers = max(1, min(workers, mp.cpu_count()))
    settings = dict(reactor=asdict(hyd.SMR_REACTOR), configuration="v3_reproduction",
                    T_grid_K=T_GRID.tolist(), flow_grid_ml_min_RTP=FLOW_GRID.tolist(),
                    T_axial_K=T_AXIAL.tolist(), flow_axial_ml_min_RTP=FLOW_AXIAL,
                    z_axial_mm=Z_MM.tolist(), workers=workers, solver=hyd.SOLVER,
                    rtol=1e-8*hyd.TOL_SCALE, atol_species=1e-16*hyd.TOL_SCALE,
                    atol_T=1e-6*hyd.TOL_SCALE, KP_MODEL=hyd.KP_MODEL, h2_seed=1e-6,
                    experimental_benchmark="853 K, 120 mL/min; 80 mL/min not validated here")
    with record_run(out, "contours", settings):
        print(f"using {workers} worker process(es)")
        with Pool(workers) as pool:
            fig_axial(pool, str(out))
            fig_operating(pool, str(out))
    print("done: contours written to", out)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    main(a.out, a.workers)
