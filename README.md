# SMR laboratory-reactor model — Costamagna et al. (2020), *Energies* 13, 2624

Python reproduction of the paper's 1-D steady-state SMR model (PBR → P-R1 → P-R2), its validation against Fig. 4
(853 K, 120 mL/min, 5 % CH₄ / 20 % H₂O / 75 % He), and extra contour analysis.

## Folder contents (keep together)

| File | Purpose |
|---|---|
| `hyd.py` | The model (kinetics, heat transfer, 3-section reactor chain) + basic plots/CSV |
| `validation.py` | Fig. 4 overlay, error table vs experiment, curve comparison, HMMC equilibrium check |
| `compare_h_length.py` | Tube- vs particle-diameter heat-transfer length scale vs the paper |
| `hyd_contours.py` | Operating-window contours + axial maps (CH₄, H₂, CO, CO₂, H₂O, ΔT) |
| `run_all.py` | Runs everything in order |
| `data/fig4_experimental_points.csv` | Experimental symbols of Fig. 4 (read from the figure; column `digitization_uncertainty`) |
| `data/fig4_digitized_curves.csv` | Paper's simulated curves of Fig. 4 (digitised): 18 curve sets = 6 variables × 3 catalysts |

## Setup (once)

```
pip install numpy scipy matplotlib
```
Python 3.9+ recommended. No other packages are needed.

## How to run

Open a terminal in the folder, then:

```
python run_all.py                  # everything (several minutes, mostly the contours), outputs in ./results
python run_all.py --workers 2      # fewer parallel processes for the contour stage (default 4)
```

or step by step:

```
python hyd.py                      # model: 3 catalysts, 853 K, 120 mL/min  -> 01_smr_results.png, 02_axial_fields_*.png, smr_summary.csv
python hyd.py --seed-study         # + H2-seed sensitivity (numerical regularisation only)
python hyd.py --validate           # + runs validation.py
python validation.py               # validation_table.csv/.md, validation_fig4_overlay.png, validation_equilibrium.csv
python compare_h_length.py         # h_length_comparison.csv/.png
python hyd_contours.py --workers 4 # contours_operating_window.png, contours_axial_<catalyst>.png (several minutes; depends on hardware)
```

Useful options of `hyd.py`: `--catalyst HMMC|57-4Q|25-4Q`, `--T 853`, `--flow 120`, `--pr2-mm 3`, `--with-pr1`,
`--h-length tube|particle`, `--out folder`. Every script accepts `--out folder` (except where noted in its header).

## What is modelled (frozen, do not retune)

* PBR 12 mm (kinetics + energy balance) → P-R1 50 mm (gas-phase WGS) → P-R2 40 mm (wall WGS, 523 K wall).
* Xu–Froment kinetics with the catalyst factors ζ₁–ζ₆ of Table 3; equilibrium constants model "D"
  (Table 4 with the E-values swapped between the columns, because the printed Table 4 is not physical).
* Heat transfer: Singhal Nu, Zehner–Schlünder k_r (corrected sign), **tube diameter** as length scale in Nu.
* A 1e-6 H₂ seed removes the p_H₂ = 0 singularity — numerical regularisation, not a physical feed component.
* Fig. 4-style plots show PBR + the first 3 mm of P-R2 only; P-R1 is simulated but not drawn.
  Outlet values (validation, operating window) use the full 40 mm P-R2, where the experimental sample is taken.

## Things to state in the report

1. **Experimental data.** The paper gives no numerical table. The experimental points were read manually from the
   Fig. 4 image. The column `digitization_uncertainty` (about ±0.0005 for CH₄/CO/CO₂, ±0.001 for H₂/H₂O, ±0.002 for the
   HMMC H₂/H₂O diamonds partly hidden by the 57-4Q circles) is *my figure-reading error*, **not** the experimental
   measurement uncertainty of Costamagna et al. (which is not published).
2. **Experimental symbols.** They are the P-R2 *outlet* compositions (40 mm). As in the paper's Fig. 4 they are drawn
   at the end of the displayed 3-mm P-R2 segment, which does not mean they were measured 3 mm into P-R2.
3. **Heat-transfer length scale.** The paper text specifies the catalyst particle diameter as the characteristic length
   in the PBR Nusselt number. Both the particle and the tube diameter were evaluated (`compare_h_length.py`). The tube
   diameter reproduces the published Fig. 4 profiles substantially better (species-curve RMSE ≈ 0.0006 vs ≈ 0.0017 (mean over catalysts);
   temperature RMSE 1.6–3.5 K vs 1.7–17 K (tube better for HMMC and 57-4Q, similar for 25-4Q); commercial-catalyst temperature dips ≈ 22 K / 15 K vs 20 K / 15 K in the text,
   particle diameter gives 13 K / 9 K) and was therefore retained for the Fig. 4 reproduction. This is a documented deviation
   from the paper's wording. Table 5 of the paper cannot be reproduced with either choice (different gas properties).
4. **Two kinds of agreement.** (A) with the paper's *published simulation* — excellent (mole-fraction RMSE ≤ 0.0015,
   temperature RMSE ≤ 3 K); (B) with the paper's *experiment* — good except 57-4Q, where the model gives CH₄ ≈ 0.0068 against
   ≈ 0.0042 measured. The paper's own simulated curve shows the same 0.0068, so this is a property of the published
   kinetics, not a coding error; do not retune the kinetics to hide it.
5. **Digitised curve coverage.** All 18 curve sets exist, but the 57-4Q line is partly hidden behind the others: 57-4Q CO₂
   has no points for z = 2.7–4 mm, 57-4Q temperature only for z = 0.15–2.4 mm (mainly the inlet dip), HMMC temperature
   z = 0.15–8.5 mm. `h_length_comparison.csv` lists the number of curves (`n_species_curves`) and temperature points
   (`n_T_points`) behind each RMSE.
6. **Contours.** The model is 1-D: the "axial contours" are z-versus-furnace-temperature maps, not (r, z) spatial contours.
   "CO / converted CH₄" is this work's own metric. Operating-window values outside 853 K and 80/120 mL/min are
   extrapolations (red stars = the validated conditions).
7. **Equilibrium.** HMMC leaves the PBR at thermodynamic equilibrium; H₂ is within 1 % of its final value from ≈ 1.6 mm and
   the maximum composition difference to equilibrium at 2 mm is ≈ 0.0005 mole fraction (paper: "after about 2 mm").

## Expected results (853 K, 120 mL/min)

HMMC exits the PBR at equilibrium (H₂ ≈ 0.160, CH₄ ≈ 0.002); 57-4Q H₂ ≈ 0.145, CH₄ ≈ 0.0068; 25-4Q H₂ ≈ 0.119,
CH₄ ≈ 0.0151. Atom balances close to ~1e-15. Temperature dips: HMMC ≈ −110 K, 57-4Q ≈ −22 K, 25-4Q ≈ −15 K.

## Troubleshooting

* `FileNotFoundError ... data/...` — keep the `data/` folder next to `validation.py` (validation is mandatory).
* `ModuleNotFoundError: hyd` — run the scripts from inside this folder.
* `FAILED operating point ...` in `hyd_contours.py` — the solver failed at that condition; the message gives the
  catalyst, temperature and flow (the error is deliberately not swallowed).
* Contour stage too slow — lower `--workers`/use a faster machine; the other stages take well under a minute.
