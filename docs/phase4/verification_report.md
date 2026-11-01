# Phase 4 — full numerical verification and final simulation results

**Status: PARTIAL.** Full workflow, reference regression, adaptive-solver convergence,
independent Radau comparison, all-catalyst seed sensitivity, state conservation,
refined energy quadrature and corrected axial contour checks pass. The preserved
400-point V3 HMMC output does **not** resolve its temperature minimum within the
predeclared 0.1 K criterion; 800 points also fail, while 1600 pass. These failures
remain recorded. This is an output-sampling limitation, not a change to governing
physics or a demonstrated failure of the adaptive integration.

## Scope, protocol and evidence integrity

Acceptance criteria were committed before execution as `357cab1`; see
[acceptance_criteria.md](acceptance_criteria.md). Primary configuration remains
`v3_reproduction`; kinetics, equilibrium interpretation, transport laws, balances,
feed, pressure and geometry are unchanged. No parameters were fitted.

Python 3.11.1, NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10.3; the isolated
`.venv-phase1` environment was reused without changing its pinned dependencies.
[Dependency check](dependency_check.txt) passes. This is verification on one
macOS/ARM environment, not a cross-platform or newly installed environment test.

All new outputs are under `phase4_runs/`. Baseline `data/` and `results/` were
compared byte-for-byte with `v3-baseline`, and Phase 1–3 outputs were preserved.
Manifests record commands, solver/settings, dependency versions, source/input
hashes, Git revision and dirty state, timestamps and generated-output hashes.
Runs occurred during development; their recorded dirty state is retained, not
rewritten to suggest a clean checkout. The final full workflow was executed at
`73c9258`, before later verification-only CLI exit-status changes.

`COMPLETED` means the calculation executed; `verification_status` is its separate
PASS/FAIL verdict. Earlier verification runners returned process code 0 after
writing a FAIL verdict. Current numerical/contour verification commands return
nonzero for failed or missing verification verdicts; automated tests exercise
this contract. Earlier manifests/logs are preserved without changing their status.

## Verification categories

| Category | Status | Evidence |
|---|---|---|
| Complete `run_all.py`, four stages | PASS | Original and final workflows completed; final wall time 132.7 s. |
| Required files and manifest hashes | PASS | 10 PNGs, 7 CSVs, 7 NPZs, 5 execution manifests, 660 contour state-check records. |
| Default regression | PASS | Four historical CSVs: zero numeric difference after documented schema mapping; all three full profile arrays exactly equal Phase 1. |
| LSODA tolerance refinement | PASS | Factors 1, 0.1, 0.01; all catalysts/sections. |
| Independent stiff solver | PASS | Radau factor 0.1 versus LSODA factor 0.01, identical equations and dense comparison grid. |
| Hydrogen seed sensitivity | PASS | 1e-6, 1e-8, 1e-10, 1e-12, all catalysts/sections. |
| Output-grid states | PASS | PBR/P-R1 400/800/1600, P-R2 2000/4000/8000; shared-point solutions stable. |
| Output-grid minimum-temperature sampling | PARTIAL | HMMC 400/800 fail 0.1 K criterion; 1600 passes. Other catalysts pass. |
| Elemental/helium and admissibility | PASS | All 24 saved sensitivity cases and 660 contour cases. |
| Independent integrated energy | PASS after refinement | Nine sections pass unchanged thresholds; initial failed quadrature evidence retained. |
| Operating-window contours | PASS at checked points | All nodes physically valid; 15 independent node checks and 3 off-grid checks pass. |
| Axial parametric maps | PASS after sampling correction | Stored 3200-sample maps checked against Radau at 700/850/1100 K for all catalysts. |
| Automated tests | PASS | 21 tests; see test_results.txt. |
| Final important figures | PASS visual review | All 10 final figures inspected; inventory below. |

## Numerical convergence and seed sensitivity

At the Figure 4 conditions (853 K, 120 mL/min), compare all species and temperature
through PBR, P-R1 and full P-R2. Thresholds: absolute species difference <=2e-5,
temperature <=0.1 K, full-outlet methane conversion <=1e-5. Relative tolerance is
1e-8 and absolute tolerances are 1e-16 mol/s and 1e-6 K before multiplication
by the listed factor. Dense verification uses a combined inlet-logarithmic and
linear grid, distinct from output sampling. Radau is an independent integration
algorithm, not an independent kinetic or thermodynamic model.

| Catalyst | Comparison | Max species difference | Max T difference (K) | Outlet conversion difference | Status |
|---|---|---:|---:|---:|---|
| HMMC | LSODA 1 vs 0.01 | 7.7188909e-10 | 2.5253677e-06 | 4.75886e-12 | PASS |
| HMMC | Radau 0.1 vs LSODA 0.01 | 1.4949618e-11 | 1.8776348e-07 | 2.8532732e-13 | PASS |
| HMMC | Seed 1e-6 vs 1e-12 | 9.99998e-07 | 0.0017499121 | 6.3816804e-07 | PASS |
| 57-4Q | LSODA 1 vs 0.01 | 7.8432276e-10 | 1.3817425e-05 | 4.9935966e-10 | PASS |
| 57-4Q | Radau 0.1 vs LSODA 0.01 | 2.0392368e-11 | 4.1061685e-07 | 2.9540814e-11 | PASS |
| 57-4Q | Seed 1e-6 vs 1e-12 | 9.99998e-07 | 0.0017499899 | 1.4215371e-06 | PASS |
| 25-4Q | LSODA 1 vs 0.01 | 6.5230184e-10 | 1.3680249e-05 | 8.5488505e-11 | PASS |
| 25-4Q | Radau 0.1 vs LSODA 0.01 | 1.681328e-11 | 4.0925102e-07 | 2.1545654e-11 | PASS |
| 25-4Q | Seed 1e-6 vs 1e-12 | 9.99998e-07 | 0.0017499375 | 2.4620705e-06 | PASS |

Intermediate tolerance/seed comparisons also pass; all values are in
[comparison_results.csv](../../phase4_runs/2026-11-01-numerics-refined/comparison_results.csv).
The approximately 1e-6 maximum seed composition difference includes the changed
inlet itself; the V3 normalization convention is unchanged.

Output-grid common-point differences are <=1.12e-16 in species, <=1.14e-13 K,
and zero in final conversion. Output counts do not constrain adaptive solver
steps. Extrema can still fall between output points:

| Catalyst | 400-point Tmin error (K) | 800-point error (K) | 1600-point error (K) | Dense Tmin (K) | Dense minimum z (mm) |
|---|---:|---:|---:|---:|---:|
| HMMC | 0.408216099 | 0.107844142 | 0.0132649272 | 743.092652272 | 0.039362929 |
| 57-4Q | 2.97730594e-06 | 6.98120402e-07 | 4.47204445e-06 | 831.424389363 | 0.240398581 |
| 25-4Q | 0.0110558523 | 0.000101809522 | 8.13734117e-05 | 837.626409157 | 0.223955260 |

HMMC's original sampled minimum is about 743.500868 K, versus dense 743.092652 K.
Use the dense result for a precise temperature-dip statement (109.907348 K from
853 K). The original summary/figures remain unchanged for regression; their
sampled extrema must not be presented as converged extrema. No default grid
or governing equation was changed. The overall numerical experiment manifest
remains **FAIL** because its HMMC 400/800 sampling checks fail.

## Physical consistency

Saved state checks test raw molar flows and temperature without clipping or
projection. Across the final 660 contour cases, maximum relative C/H/O/He residuals
are 6.835e-15 / 6.075e-15 / 6.088e-15 / 0. For the 24 sensitivity cases they are
5.126e-15 / 6.075e-15 / 6.835e-15 / 0; maximum mole-fraction normalization error
is 2.22e-16. Minimum flow is 0 mol/s; all states are finite, positive in absolute
temperature, and have admissible conversion. Final contour temperature range is
523.000312–1100.596324 K; wall cooling and exothermic WGS need not keep gas
below the furnace setpoint. Operating-window outlet conversion spans
0.09935877–0.99993174. Helium flow is exactly conserved in the saved arrays.
**Total molar flow is not required to be conserved.** Checks cover saved states,
not every internal solver stage or every continuous point in the operating domain.

## Independent integrated energy check

Coordinate V=A*z [m³]; Fi [mol/s]; hi [J/mol]; cp_i [J/mol/K]; U [W/m²/K].
Wall source is U*(4/d)*(Twall-T) [W/m³], positive into the gas. Reaction source is
`-sum(hi*dFi/dV)` [W/m³], negative for net endothermic reforming. Independently
check endpoint total enthalpy change `Delta(sum Fi*hi)=Qwall`, and temperature
line integral `integral(sum Fi*cp_i)dT=Qwall+Qreaction`, all in W. Explicit
stoichiometric source expressions are evaluated separately from the production
stoichiometric matrix; enthalpy and heat capacities still use the same property
model. This is an integration/sign/unit consistency check, not an independent
validation of thermodynamic correlations.

Use LSODA tolerance factor 0.001 for this check, independent quadrature with
15999 combined logarithmic/linear nodes and a refinement comparison. Relative
residual denominator is max(|Qwall|,|Qreaction|,|Delta H|,1e-6 W); threshold 1e-4,
quadrature-change threshold 1e-5. Both were fixed before execution.

| Catalyst | Section | Wall heat (W) | Reaction heat (W) | Total-enthalpy residual (W) | Sensible-balance residual (W) | Max relative residual | Status |
|---|---|---:|---:|---:|---:|---:|---|
| HMMC | PBR | 0.849609994 | -0.842089082 | -6.79821e-11 | -1.40983e-09 | 1.65938e-09 | PASS |
| HMMC | PR1 | -1.44290232e-12 | -2.46935553e-23 | 1.57213e-12 | 1.57209e-12 | 1.57213e-06 | PASS |
| HMMC | PR2 | -0.817144688 | 0.00885124514 | 1.22941e-10 | 1.97446e-10 | 2.41629e-10 | PASS |
| 57-4Q | PBR | 0.753657313 | -0.756806039 | 3.40927e-11 | -2.64507e-11 | 4.50482e-11 | PASS |
| 57-4Q | PR1 | 0.00390341037 | 5.37474684e-12 | -2.80899e-11 | -2.80893e-11 | 7.19624e-09 | PASS |
| 57-4Q | PR2 | -0.813743338 | 0.00989698533 | 1.22833e-10 | 2.71842e-10 | 3.34064e-10 | PASS |
| 25-4Q | PBR | 0.590315855 | -0.595001458 | 1.68153e-11 | -1.69024e-11 | 2.84073e-11 | PASS |
| 25-4Q | PR1 | 0.00511912684 | 5.43376311e-12 | -4.08721e-11 | -4.08725e-11 | 7.98427e-09 | PASS |
| 25-4Q | PR2 | -0.804900769 | 0.00879119598 | 1.18624e-10 | 3.10317e-10 | 3.85534e-10 | PASS |

Largest absolute residual is 1.410e-9 W; largest normalized residual is
1.573e-6 in the near-zero-duty HMMC P-R1, using the predeclared 1e-6 W floor.
Maximum quadrature refinement change is 4.672e-9 relative.

The first energy attempt failed HMMC PBR (about 0.003018 W sensible residual)
and its near-zero-duty P-R1. The rate denominator at seeded inlet produces an
extremely narrow initial reaction/heat layer. Extending log quadrature to L*1e-18
from L*1e-12 resolves it, and tighter integration resolves the near-zero P-R1
cancellation. Thresholds, source physics and seed were not relaxed or retuned.
Initial failed results remain in `2026-11-01-numerics/energy_balance.csv`.

## Full workflow and contour verification

The first full workflow completed at `phase4_runs/2026-11-01-workflow/`.
Independent axial-map checks then exposed interpolation errors from 400 source
samples: up to 1.4943 K and 7.019e-4 in mole fraction in the checked cases.
Testing 1600 samples still left high-temperature composition failures; 3200
passed all nine checked cases. Only `hyd_contours._axial` now requests 3200
source samples. The adaptive solver, its tolerances and `hyd.run` default remain
unchanged. The final full workflow was rerun after this numerical plotting fix
at `phase4_runs/2026-11-01-final-workflow/`; all four child processes exit 0.

Operating grid: catalysts × 17 furnace temperatures (700–1100 K) × 11 feed rates
(60–200 mL/min at RTP) × 3 dimensionless fields = (3,17,11,3). Axial grid per
catalyst: 33 furnace temperatures × 6 fields × 241 positions (0–12 mm) =
(33,6,241). Species are mole fractions; the sixth axial field is T-Tf [K].
The 561 operating and 99 axial cases all generated passing raw-state records.
No solver failures/invalid regions occurred. CO/converted-CH4 is an explicitly
labelled derived metric; its below-5%-conversion mask is saved and has zero
masked points in this domain. Physical checks precede plotting; invalid states
raise instead of becoming hidden NaNs. Raw NPZ arrays are saved before plots.

Operating checks independently solve the four corners and central node for each
catalyst with 100x tighter tolerances (15 node checks): maximum field difference
3.023e-8, below 2e-5. At off-grid 853 K/120 mL/min, bilinear interpolation versus
a direct solve differs by at most 0.000894318 / 0.001975059 / 0.000207186 for
HMMC/57-4Q/25-4Q, below the separately declared 0.005 interpolation limit.
The benchmark is not an exact node of the operating sweep.

Axial reconstruction agrees with the matching default interpolant to roundoff.
The stronger check evaluates Radau directly at stored z positions at 700, 850,
and 1100 K for every catalyst:

| Catalyst | Corrected map vs Radau: max species error | Max T error (K) | Max old-to-new species change across all maps | Max old-to-new T change (K) |
|---|---:|---:|---:|---:|
| HMMC | 9.40319226e-06 | 0.00666012348 | 0.00069248186 | 0.619508082 |
| 57-4Q | 5.95252732e-06 | 0.0195734353 | 0.00045328623 | 1.4746774 |
| 25-4Q | 5.03141359e-06 | 0.0165366341 | 0.000381833128 | 1.24032603 |

Operating raw arrays and all reference predictions are unchanged exactly.
Map corrections affect interpolation only and improve numerical representation;
they do not introduce physical assumptions. The 0.05 mm display grid still cannot
locate sub-cell thermal extrema: use the dense numerical minima above. Nine
independent axial cases do not establish an error bound over all 33 temperatures
or between plotted nodes. Axial maps are position-versus-operating-temperature
parametric maps, not radial CFD contours. Sweeps are predictions beyond the
fixed validated benchmark, including 80 mL/min.

## Recalculated Figure 4 and experimental metrics

Compare digitized PBR numerical curves only at supplied points: species z>0.05 mm,
temperature z>0.15 mm; visibility/coverage gaps remain. Mean species RMSE is the
unweighted mean of five species RMSEs, not a pooled error. The historical
heat-length script additionally uses z<5 mm for temperature, a different domain.
Experimental metrics use CH4/H2/CO/CO2/H2O at the **full 40 mm P-R2 outlet**, at
853 K and 120 mL/min. MAPE is mean absolute relative error ×100 across those five
figure-read compositions. Digitization uncertainty is not measurement uncertainty.

| Catalyst | Mean species curve RMSE | Temperature RMSE, all retained points (K) | Experimental MAE | Experimental RMSE | Experimental MAPE |
|---|---:|---:|---:|---:|---:|
| HMMC | 0.0006296162264 | 1.2354092 | 0.001217200332 | 0.00157951952 | 6.3406102% |
| 57-4Q | 0.0005633421608 | 2.88919175 | 0.00394211702 | 0.005094651504 | 16.575256% |
| 25-4Q | 0.0005136485294 | 2.97130116 | 0.001263318078 | 0.001555044141 | 3.8500005% |

The inlet-domain temperature RMSEs remain 1.638515 / 2.889192 / 3.541152 K.
All historical numeric values in smr_summary.csv, validation_table.csv,
validation_equilibrium.csv and h_length_comparison.csv match V3 exactly after
the Phase 3 schema/label mapping. All saved species, temperature and flow arrays
match Phase 1 exactly. Agreement with digitized simulations and measured outlets
is reported separately. Lower RMSE does not prove physical superiority of V3.

## Figure inventory and visual inspection

All paths below are relative to `phase4_runs/2026-11-01-final-workflow/`.
All ten PNGs were opened and inspected for axes, units, legends, annotations,
visible gradients and layout. None showed a new material clipping/legibility
problem. Machine inventory includes dimensions, byte lengths and SHA-256 hashes
in `2026-11-01-final-artifact-audit-02/file_inventory.json`. Visual review is not
an additional measurement of numerical accuracy.

| Figure | Purpose |
|---|---|
| `model/01_smr_results.png` | Six reference species/temperature profiles for three catalysts. |
| `model/02_axial_fields_HMMC.png` | 1-D axial field strips, HMMC. |
| `model/02_axial_fields_57-4Q.png` | 1-D axial field strips, 57-4Q. |
| `model/02_axial_fields_25-4Q.png` | 1-D axial field strips, 25-4Q. |
| `fixed_reference/validation_fig4_overlay.png` | Published simulated profiles plus full-outlet experimental markers. |
| `h_length_sensitivity/h_length_comparison.png` | Historical two-case Nu-to-h sensitivity. |
| `contours/contours_operating_window.png` | Predicted operating window; fixed benchmark distinguished. |
| `contours/contours_axial_HMMC.png` | Refined 1-D axial-position/temperature parametric maps. |
| `contours/contours_axial_57-4Q.png` | Same for 57-4Q. |
| `contours/contours_axial_25-4Q.png` | Same for 25-4Q. |

## Execution and reproducibility

The original and corrected full workflows were each run once; the second was
justified by the confirmed contour interpolation error. Stages were not rerun
individually for an additional duplicate pipeline. Dedicated verification tools
perform additional solver calls for convergence and independent comparisons.
Exact historical invocations are recorded in every run manifest and stdout/stderr
logs are copied into [execution_logs/](execution_logs/). Use new destinations;
existing nonempty outputs are refused.

From the project root, with the pinned environment and writable Matplotlib cache:

```sh
export MPLCONFIGDIR=/tmp/hydrogen-phase4-mpl
export PYTHONDONTWRITEBYTECODE=1
.venv-phase1/bin/python run_all.py --workers 2 --out phase4_runs/replay-final-workflow
.venv-phase1/bin/python tools/verify_numerics.py --out phase4_runs/replay-numerics
.venv-phase1/bin/python tools/verify_contours.py --workflow phase4_runs/replay-final-workflow --out phase4_runs/replay-contours
.venv-phase1/bin/python tools/check_axial_resolution.py --workflow phase4_runs/replay-final-workflow --out phase4_runs/replay-axial
.venv-phase1/bin/python tools/audit_final_artifacts.py --workflow phase4_runs/replay-final-workflow --out phase4_runs/replay-artifacts
.venv-phase1/bin/python -m unittest discover -s tests -v
```

`verify_numerics.py` is expected to return **1** with the current preserved default
HMMC sampling failures; do not ignore its CSV/manifest verdict. The refined
energy experiment was executed before this exit-code-only fix and returned 0
while its manifest correctly recorded FAIL. The complete expensive experiment
was not repeated just to re-test exit handling; automated tests verify the new
verdict-to-exit-code mapping. Final contour and axial checkers were rerun after
the exit-status correction into `*-02` evidence directories. No physical or
reference-prediction changes occurred between those revisions.

The artifact audit assumes preserved Phase 1 data, the annotated V3 Git tag and
the original Phase 4 workflow are present for comparison. It is a repository
verification command, not a standalone model execution entry point.

## Changes and commits

| Commit | Scope |
|---|---|
| `357cab1` | Predeclared acceptance criteria. |
| `f50b4a4` | Raw reference profiles, contour grids and per-case physical evidence. |
| `bde4306` | Independent numerical/energy verification tools. |
| `2e8782c` | Refined inlet quadrature and explicitly measured output-grid differences. |
| `73c9258` | Contour source-sampling correction (400 to 3200), focused tests. |
| `8d082bb` | Distinguish stored resolution from candidate resolutions. |
| `958b893` | Numerical verification failures propagate to CLI status. |
| `587198b` | Final artifact audit, all checker exit contracts and additional tests. |

Production changes: hyd.py saves profiles; hyd_contours.py records raw/evidence
outputs, rejects invalid states and refines axial resampling. numerical_checks.py
adds checks without altering state. New tools are verify_numerics.py,
verify_contours.py, check_axial_resolution.py and audit_final_artifacts.py;
focused tests are in tests/test_numerical_checks.py. README and this report link
current figures/evidence. No reference data or baseline outputs were modified.
The final evidence/documentation commit is identified in Git history and the
phase-completion reply; it cannot embed its own final hash.

## Remaining limitations and next phase

- The preserved default HMMC output grid misses its narrow minimum. Use dense
  evidence for quantitative minima; approving a future default sampling change
  is separate from preserving the V3 regression record.
- Table 5 versus printed transport equations/properties remains scientifically
  ambiguous; the reproduction configuration is not physically proven superior.
- Existing equilibrium interpretation, neglect of dispersion/pressure drop/radial
  fields, property extrapolation and catalyst-model applicability remain model
  limitations. Numerical convergence does not remove them.
- Energy checks use the model's own property functions; equilibrium checks share
  equilibrium constants with kinetics. Neither is independent physical validation.
- Figure-read data and incomplete digitization coverage constrain validation;
  no supplied 80 mL/min experimental benchmark is validated here.
- Contour checks are finite-domain sampling, not a rigorous continuous error bound.
- Matplotlib/Pyparsing deprecation warnings remain; no dependency upgrades were
  introduced. The test log's synthetic invalid-state failure message is deliberate
  and its rejection test passes; it is not a hidden production solver failure.
- The actual university assignment PDF has not been checked in this phase.
  The academic report, presentation, compliance checklist and submission package
  belong to Phase 5 and have not been started. No merge into main is authorized.

**Stop after Phase 4. Request approval before Phase 5.**
