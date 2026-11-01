"""
validation.py - validation of the SMR model (hyd.py) against Costamagna et al., Energies 2020, 13, 2624, Fig. 4.

CASE   5 % CH4 / 20 % H2O / 75 % He, 120 mL/min (RTP), furnace 853 K (GHSV 2.1e4 1/h), 1 atm.

WHAT IT DOES
  1. runs the model for HMMC, 57-4Q and 25-4Q;
  2. compares the outlet composition with the EXPERIMENTAL symbols of Fig. 4 (table, abs. and % error,
     MAE / RMSE / mean % error);
  3. compares the simulated curves with the PAPER'S SIMULATION curves of Fig. 4 (RMSE per species);
  4. draws the Fig. 4 reproduction: own simulation (dashed) + paper simulation (thick) + experiments (symbols);
  5. checks that HMMC reaches thermodynamic equilibrium inside the PBR (as stated in the paper).

DATA (mandatory, in ./data/ next to this file; the script stops with an error if they are missing)
  fig4_experimental_points.csv  experimental symbols, READ FROM THE FIGURE because the paper gives no table
                                (5x zoom, calibrated axes). Column digitization_uncertainty is MY estimate of the
                                figure-reading error (+-0.0005 CH4/CO/CO2, +-0.001 H2/H2O, +-0.002 for the HMMC
                                H2/H2O diamonds that are partly hidden by the 57-4Q circles). It is NOT the
                                experimental measurement uncertainty of Costamagna et al., which is not published.
  fig4_digitized_curves.csv     the paper's simulated curves (18 curve sets = 6 variables x 3 catalysts), digitised by
                                colour from the PDF figure. Coverage is limited to visible PBR portions; some lines are hidden behind
                                another: 57-4Q CO2 has a gap between z = 2.5512-4.0866 mm, 57-4Q T spans z = 0.148-1.7946 mm
                                (the black line lies under the green/red ones), HMMC T for z = 0.15-8.5 mm.

TWO KINDS OF AGREEMENT ARE REPORTED SEPARATELY
  (A) model vs the paper's PUBLISHED SIMULATION (does this code reproduce the published model?)
  (B) model vs the paper's EXPERIMENT (does the model describe the laboratory data?)
  The published model itself deviates from experiment (notably 57-4Q CH4), so (A) can be excellent while (B) is not.

OUTPUT   validation_table.csv, validation_table.md, validation_equilibrium.csv, validation_fig4_overlay.png
RUN      python validation.py          (or: python hyd.py --validate)
"""
import csv
import os
import numpy as np
import hyd as h

HERE = os.path.dirname(os.path.abspath(__file__))
SPECIES = ["CH4", "H2", "CO", "CO2", "H2O"]                       # species with experimental points
IDX = {"CH4": 0, "H2O": 1, "CO": 2, "CO2": 3, "H2": 4}           # row of each species in the model state
CATS = ["HMMC", "57-4Q", "25-4Q"]
TF, FLOW = 853.0, 120.0                                           # the Fig. 4 operating condition


# --------------------------------------------------------------------------------------------------
# Data loading (mandatory: raises FileNotFoundError instead of silently skipping validation)
# --------------------------------------------------------------------------------------------------
def _need(fname):
    path = os.path.join(HERE, "data", fname)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"validation needs {path} - keep the data/ folder next to validation.py")
    return path


def load_experiment():
    """returns {catalyst: {species: (value, uncertainty)}} from fig4_experimental_points.csv"""
    out = {}
    with open(_need("fig4_experimental_points.csv"), newline="") as stream:
        for r in csv.DictReader(stream):
            out.setdefault(r["catalyst"], {})[r["species"]] = (float(r["mole_fraction"]), float(r["digitization_uncertainty"]))
    return out


def load_paper_curves():
    """returns {(species, catalyst): array[[z_mm, value], ...]} - paper simulation curves in the PBR (T in K)"""
    out = {}
    with open(_need("fig4_digitized_curves.csv"), newline="") as stream:
        for r in csv.DictReader(stream):
            out.setdefault((r["species"], r["catalyst"]), []).append((float(r["z_mm"]), float(r["value"])))
    return {k: np.array(v) for k, v in out.items()}


def exit_composition(res):
    """mole fractions at the physical P-R2 outlet (40 mm): this is where the experimental sample is taken"""
    y = h.mole_fractions(res["PR2"][1])[:, -1]
    return {s: float(y[i]) for s, i in IDX.items()}


# --------------------------------------------------------------------------------------------------
# Main validation routine
# --------------------------------------------------------------------------------------------------
def _main(out):
    os.makedirs(out, exist_ok=True)
    P = lambda f: os.path.join(out, f)
    EXP, CURVES = load_experiment(), load_paper_curves()

    # ---- 1) run the model for the three catalysts -----------------------------------------------
    results = {c: h.run(catalyst=c, Tf=TF, flow=FLOW) for c in CATS}

    # ---- 2) outlet composition vs experiment: table + error statistics --------------------------
    rows = []
    for c in CATS:
        ym = exit_composition(results[c])
        for s in SPECIES:
            e, tol = EXP[c][s]                                     # experiment, figure-reading uncertainty
            m = ym[s]
            rows.append(dict(catalyst=c, species=s, experiment=e, model=round(m, 5), signed_error=round(m - e, 5), absolute_error=round(abs(m - e), 5),
                             rel_err_pct=round(100 * (m - e) / e, 1), digitization_uncertainty=tol,
                             within_digitization_uncertainty=bool(abs(m - e) <= tol)))
    with open(P("validation_table.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(P("validation_table.md"), "w") as f:             # same table, ready to paste into the report
        f.write("Experimental values were read from Fig. 4 (no numerical table in the paper); model = P-R2 outlet.\n\n")
        f.write("| Catalyst | Species | Experiment (Fig. 4) | This work | Signed error | Signed % error |\n|---|---|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(f"| {r['catalyst']} | {r['species']} | {r['experiment']:.4f} | {r['model']:.4f} | "
                    f"{r['signed_error']:+.4f} | {r['rel_err_pct']:+.1f} |\n")

    print("(B) MODEL vs EXPERIMENT at the P-R2 outlet (experiment values read from Fig. 4)")
    print(f"{'catalyst':9s}{'species':8s}{'exp':>9s}{'model':>9s}{'signed err':>10s}{'% err':>8s}   within digitization uncertainty")
    for r in rows:
        print(f"{r['catalyst']:9s}{r['species']:8s}{r['experiment']:9.4f}{r['model']:9.4f}"
              f"{r['signed_error']:+10.4f}{r['rel_err_pct']:+8.1f}   {r['within_digitization_uncertainty']}")
    print("\nError statistics of the outlet composition (MAE, RMSE, mean |% error|):")
    for c in CATS:
        e = np.array([r["signed_error"] for r in rows if r["catalyst"] == c])
        p = np.array([abs(r["rel_err_pct"]) for r in rows if r["catalyst"] == c])
        print(f"  {c:6s} MAE={np.mean(abs(e)):.4f}  RMSE={np.sqrt(np.mean(e**2)):.4f}  mean|%err|={np.mean(p):.1f}")

    with open(P("experimental_metrics.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["catalyst", "MAE", "RMSE", "MAPE_pct", "n_species"])
        writer.writeheader()
        for c in CATS:
            ym = exit_composition(results[c])
            errors = np.array([ym[s]-EXP[c][s][0] for s in SPECIES])
            percentages = [100*abs(ym[s]-EXP[c][s][0])/EXP[c][s][0] for s in SPECIES]
            writer.writerow(dict(catalyst=c, MAE=np.mean(abs(errors)), RMSE=np.sqrt(np.mean(errors**2)),
                                 MAPE_pct=np.mean(percentages), n_species=len(SPECIES)))
    with open(P("digitization_coverage.csv"), "w", newline="") as f:
        writer = csv.writer(f); writer.writerow(["catalyst", "variable", "n_points", "z_min_mm", "z_max_mm", "largest_gap_mm"])
        for (s, c), d in sorted(CURVES.items()):
            zz = np.sort(d[:,0])
            writer.writerow([c, s, len(d), zz.min(), zz.max(), np.diff(zz).max() if len(zz)>1 else 0])
    curve_metrics = []

    # ---- 3) curves vs the paper's own simulation (PBR, z > 0.05 mm) ------------------------------
    print("\n(A) MODEL vs the paper's PUBLISHED SIMULATION: RMSE of the PBR curves (mole fraction; T in K).")
    print("    (n = number of digitised curve sets used per catalyst, of 6 possible)")
    for c in CATS:
        z, y, T, b = h.profiles(results[c])
        n = z <= b["PBR_end"]
        parts = []; ncurves = 0
        for s in SPECIES + ["T"]:
            if (s, c) not in CURVES:
                continue
            ncurves += 1
            d = CURVES[(s, c)]; d = d[d[:, 0] > 0.15 if s == "T" else d[:, 0] > 0.05]
            mine = np.interp(d[:, 0], z[n], (T if s == "T" else y[IDX[s]])[n])
            rmse = np.sqrt(np.mean((mine - d[:, 1]) ** 2))
            parts.append(f"{s}={rmse:.4g}")
            curve_metrics.append(dict(catalyst=c, variable=s, n_points=len(d), z_min_mm=d[:,0].min(),
                                      z_max_mm=d[:,0].max(), RMSE=rmse, units="K" if s=="T" else "mole fraction"))
        print(f"  {c:6s} (n={ncurves}) " + "  ".join(parts))

    with open(P("published_curve_metrics.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(curve_metrics[0]))
        writer.writeheader(); writer.writerows(curve_metrics)

    # ---- 4) HMMC equilibrium check ----------------------------------------------------------------
    eq = h.equilibrium([0.05, 0.20, 0.0, 0.0, 0.0, 0.75], TF)    # feed composition, 853 K, 1 atm
    z, y, T, b = h.profiles(results["HMMC"]); n = z <= b["PBR_end"]
    print("\nHMMC PBR outlet vs thermodynamic equilibrium at 853 K (paper: equilibrium after about 2 mm):")
    with open(P("validation_equilibrium.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["species", "HMMC_PBR_exit", "equilibrium", "signed_difference"])
        for s in SPECIES:
            a_, b_ = y[IDX[s]][n][-1], eq[IDX[s]]
            w.writerow([s, round(a_, 5), round(b_, 5), round(a_ - b_, 5)])
            print(f"  {s:4s} model={a_:.4f}  equilibrium={b_:.4f}  diff={a_ - b_:+.4f}")
    k = np.argmax(np.abs(y[IDX['H2']][n] - y[IDX['H2']][n][-1]) < 0.01 * y[IDX['H2']][n][-1])
    print(f"  H2 within 1 % of its PBR-exit value from z = {z[n][k]:.2f} mm onward")
    # composition at 2 mm (the paper: "thermodynamic equilibrium after about 2 mm") vs equilibrium
    d2 = max(abs(np.interp(2.0, z[n], y[IDX[s]][n]) - eq[IDX[s]]) for s in SPECIES)
    print(f"  at z = 2 mm: max |composition - equilibrium| over CH4,H2,CO,CO2,H2O = {d2:.5f} mole fraction")
    with open(P("validation_equilibrium.csv"), "a", newline="") as f:
        csv.writer(f).writerow(["max_abs_diff_at_2mm", round(d2, 6), "", ""])

    # ---- 5) Fig. 4 reproduction plot ----------------------------------------------------------------
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"HMMC": "tab:red", "57-4Q": "k", "25-4Q": "tab:green"}
    mk = {"HMMC": "D", "57-4Q": "o", "25-4Q": "s"}               # marker shapes as in the paper
    fig, ax = plt.subplots(2, 3, figsize=(15, 8)); ax = ax.ravel()
    for a, s in zip(ax, ["CH4", "H2", "CO", "CO2", "H2O", "T"]):
        for c in CATS:
            z, y, T, b = h.profiles(results[c])                     # PBR + first 3 mm of P-R2 (P-R1 skipped)
            if (s, c) in CURVES:                                    # thick = paper simulation (PBR part)
                d = CURVES[(s, c)]; a.plot(d[:, 0], d[:, 1], color=col[c], alpha=0.25, lw=6)
            a.plot(z, T if s == "T" else y[IDX[s]], "--", color=col[c], lw=1.5)    # dashed = this work
            if s != "T":                                            # symbol = experiment at the sampling point
                a.plot([z[-1]], [EXP[c][s][0]], mk[c], mfc="none", mec=col[c], ms=9, mew=1.8)
        a.axvline(b["PBR_end"], ls="--", c="k", lw=1)               # PBR / P-R2 boundary
        top = 1.01
        a.text(b["PBR_end"] / 2, top, "PBR", ha="center", va="bottom", fontsize=9, transform=a.get_xaxis_transform())
        a.text(b["PBR_end"] + 1.5, top, "P-R2", ha="center", va="bottom", fontsize=9, transform=a.get_xaxis_transform())
        a.set_ylabel("temperature [K]" if s == "T" else "mole fraction"); a.set_title(s if s != "T" else "Temperature [K]", pad=25); a.set_xlabel("reactor length [mm]  (PBR | first 3 mm of P-R2)")
    # one shared legend below the panels (thick = paper simulation, dashed = this work, symbol = experiment)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color="k", lw=6, alpha=0.25, label="paper simulation (digitised)"),
                        Line2D([], [], color="k", ls="--", label="this work"),
                        *[Line2D([], [], color=col[c], marker=mk[c], mfc="none", ls="", label=f"{c}: experiment (digitized)") for c in CATS]],
               loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.035), fontsize=9)
    fig.suptitle("Fig. 4 reproduction: 853 K, 120 mL/min (red HMMC, black 57-4Q, green 25-4Q)")
    fig.text(0.5, 0.012, "Experimental symbols are the P-R2 OUTLET compositions (40 mm); as in the paper's Fig. 4 they are drawn at the end of "
             "the displayed 3-mm P-R2 segment.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.12, 1, .97)); fig.savefig(P("validation_fig4_overlay.png"), dpi=130, bbox_inches="tight", pad_inches=.2); plt.close(fig)
    print(f"\nWrote validation_table.csv/.md, validation_equilibrium.csv, validation_fig4_overlay.png to {os.path.abspath(out)}")
    return rows


def main(out=None):
    """Standalone fixed Figure 4 benchmark, independent of custom simulations."""
    from run_metadata import prepare_output, record_run
    from dataclasses import asdict
    out = prepare_output(out, "validation")
    settings = dict(benchmark="Figure 4", T_K=TF, flow_ml_min_RTP=FLOW,
                    reactor=asdict(h.SMR_REACTOR), pressure_bar=h.P_BAR,
                    solver=h.SOLVER, rtol=1e-8*h.TOL_SCALE, atol_species=1e-16*h.TOL_SCALE,
                    atol_T=1e-6*h.TOL_SCALE, KP_MODEL=h.KP_MODEL, h2_seed=1e-6,
                    experimental_location="P-R2 full 40 mm outlet",
                    species_curve_domain="PBR z>0.05 mm", temperature_curve_domain="PBR z>0.15 mm")
    with record_run(out, "validation", settings):
        return _main(str(out))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=None)
    main(ap.parse_args().out)
