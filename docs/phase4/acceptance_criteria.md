# Phase 4 acceptance criteria — fixed before execution

Primary configuration v3_reproduction, unchanged kinetics/equilibrium/transport,
geometry, pressure and feed. Preserve all historical evidence. Passing numerical
checks is not independent confirmation of the paper's transport interpretation.

## Execution, artifacts and regression

Run run_all.py once with 2 workers into a fresh versioned directory; each of four
stages and all manifests must complete. Require 10 figures, the four existing
numeric CSVs plus coverage/metrics CSVs, three reference profiles, raw operating
and axial contour grids, per-point physical checks and source/output hashes.
Default reference numeric CSV fields must equal V3 (zero after documented schema
mapping); full saved reference arrays must equal Phase 1. Byte equality of renamed
CSV headers/labels or newly rendered PNGs is not required.

## Convergence and sensitivities — all three catalysts at Figure 4 conditions

Tolerance factors: LSODA default, 0.1 and 0.01 times rtol/atol together. Compare
with the tightest solution on identical grids: max species mole-fraction difference
<=2e-5, max temperature difference <=0.1 K, outlet conversion difference <=1e-5.
Output grids: PBR/P-R1 400,800,1600; P-R2 2000,4000,8000. State values at common
points must satisfy the same limits; additionally sampled PBR Tmin must change
<=0.1 K against a dense-output, log-inlet/linear verification grid. Report failures
of this extrema criterion separately from adaptive-solver solution convergence.
Output sampling is not a solver step-size control.

Independent stiff solver: Radau at tolerance factor 0.1 versus LSODA factor 0.01,
same dense verification grid, same <=2e-5 species, <=0.1 K temperature,
<=1e-5 conversion limits. Solver failure is FAIL, not a missing point to ignore.
Hydrogen seed: 1e-6,1e-8,1e-10,1e-12, fixed LSODA factor 0.1; compare to smallest
seed over dense profiles and outlet: <=2e-5 absolute species mole fraction,
<=0.1 K temperature and <=1e-5 conversion. The feed is renormalized exactly as V3;
this is regularization sensitivity, not tuning the physical feed.

## Physical consistency and independent integrated energy check

Across every saved state in all sections: finite values, T>0 K, flows >=-1e-14
mol/s; mole-fraction sums within 1e-12; relative C/H/O/He residuals <=1e-8;
CH4 conversion within [-1e-8,1+1e-8]. Do not require total molar-flow conservation.

On tight dense-output solutions independently integrate wall heat:
Delta[sum_i Fi hi(T)] = integral U*(4/d)*(Twall-T) dV [W]. Also check
integral sum_i Fi cpi(T)*dT/dV dV = Qwall + Qreaction [W], where
Qreaction = -integral sum_i hi(T)*dFi/dV dV (positive for net heat release).
Use thermochemical enthalpies independently from stored temperature/composition,
and quadrature separate from the ODE stepping. Normalize residuals by
max(abs(Qwall),abs(Qreaction),abs(Delta enthalpy),1e-6 W); require <=1e-4.
Require quadrature refinement error <=1e-5 of the same scale. Preserve failed
attempts; if uniform quadrature undersamples the singular inlet, refine the
quadrature without relaxing the thresholds or altering the physics.

## Contours

Operating grid (3,17,11,3); axial grid per catalyst (33,6,241). Save coordinates,
units and un-clipped raw fields. Require 561 operating +99 axial physical-check
records, all valid; no solver errors silently discarded. NaN permitted only in
CO/converted-CH4 where conversion<0.05, with the mask explicitly recorded.
Verify all four operating corners, a central grid point and 853K/120 reference
by independently evaluating each catalyst with tighter tolerances; compare
direct node values <=2e-5, interpolation estimates at off-grid 853K/120 <=0.005
absolute (interpolation accuracy distinct from solver accuracy). Axial values
must equal independently sampled default solution interpolants at stored nodes
within 2e-5; also report difference from tighter directly evaluated profiles,
including sharp inlet effects rather than hiding them.

Inspect every important new figure. Parametric axial maps are not radial CFD.
Any invalid state, failed check or unsupported domain must be explicit in the
report. Full pipeline completion, numerical PASS and physical-model limitations
must be reported separately. New physical assumptions require user approval.
