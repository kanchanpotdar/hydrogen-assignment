# Phase 2 comparison protocol (before numerical evaluation)

Default remains V3. Only PBR transport definitions vary; kinetics, equilibrium,
solid-conductivity fit, reactor geometry, feed, solver and balances remain fixed.

Figure 4: 853 K, 120 mL/min RTP, 5/20/75 CH4/H2O/He, all three catalysts.
Use exactly the same digitized points for all configurations. Species z>0.05 mm;
temperature z>0.15 mm, with both all-available and inlet z<5 mm metrics reported.
Report per-species RMSE and their unweighted mean (not pooled RMSE), PBR thermal
dip, conversion, outlet composition and differences from V3. Experimental
MAE/RMSE/MAPE use the five measured species at the full 40 mm P-R2 outlet.
Digitization uncertainty is not experimental measurement uncertainty.

Table 5: independently compare its two SMR columns (773/1173 K). Export the
published rows. First use their tabulated rho/us/mu/cp/kf/ks directly, then use
V3 property correlations with the stated SMR feed and GHSV=2.4e4/h. Infer RTP
flow from GHSV*PBR volume/60 and disclose differences from tabulated velocity.
Fixed wall conductivity 1.6 W/m/K; ks fit retained for correlated calculations;
tabulated ks used only in the tabulated-property calculation. Never substitute
Table 5 conditions for the Figure 4 case or fit missing inputs to U.

Verification: unit-test independently expressed Re, Pr, Nu, h and U for every
configuration using fixed SI inputs and dimensional scaling. Legacy default and
particle-only interfaces must agree with named equivalents. Reject unknown
configurations and contradictory explicit CLI options. V3 summary and four
baseline CSVs must match exact saved bytes under unchanged environment. Full
precision V3 profiles must agree with Phase 1 (rtol=1e-10, atol=1e-16 mol/s for
flows, atol=1e-7 K for temperature). Alternatives: finite saved states,
minimum raw flow >=-1e-14 mol/s; relative elemental/He residual <=1e-8,
normalization error <=1e-12. All outcomes, including failures, are persisted.
No minimum fitting accuracy is imposed on alternatives. These comparisons do
not establish solver convergence or physical correctness. No default selection
based on RMSE. No new constant-solid-conductivity case in Phase 2.
