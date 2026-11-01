# Phase 1 acceptance criteria — specified before execution

Scope: unchanged V3 at 853 K, 120 mL/min RTP, all three catalysts; existing
transport comparison includes its original tube and h-only particle modes.
No alternative physics configurations, contour rerun, or full Phase 4 verification.

1. Every original Python, data, and results blob must match v3-baseline exactly.
2. hyd.py, validation.py and compare_h_length.py must exit zero in an isolated
   environment. Save commands, logs, package versions, input hashes and run status.
3. Compare CSV schema, row count and string identifiers exactly. Compare finite
   numbers with rtol=0 and the following absolute tolerances (reflecting saved precision):
   smr_summary.csv: mole fractions 5e-6, T_K 0.02 K, X_CH4 2e-4;
   validation_table.csv: mole fractions/errors 2e-5, rel_err_pct 0.2 percentage points,
   digitization_uncertainty exact and uncertainty flags exact;
   validation_equilibrium.csv: compositions/residuals 2e-5;
   h_length_comparison.csv: U 2 W/m2/K, dip/T-RMSE 0.2 K,
   species RMSE and outlet MAE 2e-5; integer counts exact.
   Empty fields and boolean labels must match. NaN/Inf are failures.
4. Report byte equality separately: any numerical-tolerance PASS is not byte equality.
5. Recompute outlet MAE/RMSE/MAPE across CH4,H2,CO,CO2,H2O. Compare curves
   only at supplied digitized points, excluding z<=0.05 mm for species and
   z<=0.15 mm for temperature. Preserve the comparison script's additional
   temperature restriction z<5 mm; report that domain separately.
6. Baseline sanity checks for all catalysts: finite states; minimum raw flow
   >=-1e-14 mol/s; mole-fraction sum error <=1e-12; relative C/H/O and helium
   flow errors <=1e-8 across all saved axial states. No total-mole conservation requirement.
7. New PNGs must exist and be visually inspected. Byte-identical PNGs are not
   required because dependency/rendering differences can affect serialization.

Reproduction is distinct from physical validation. Matching V3 does not resolve
transport ambiguity or prove convergence. Thresholds are fixed before viewing
new numerical results; failures must be retained and explained.
