# Phase 1 — Baseline reproduction and scientific specification

Status: **PASS for the defined Phase 1 scope**, 1 November 2026.
This does not certify convergence, the complete contour pipeline, or physical
correctness of every transport interpretation.

## Preservation and changes

Baseline tag: `v3-baseline`, commit `a5782cc54e89dc7117b82b076279a75452e85903`.
Execution commit: `f8e64078b53a7acf150dd3b0ebb553021fd6ecb3`.
Acceptance criteria committed before execution at `4deab9b`; harness committed
at `6597591`. The README plot addition preceded verification at `203610f`.

All five original Python scripts, both data CSVs, and every original results
file were hash-checked against the baseline before and after execution: unchanged.
No governing equations, kinetic factors, equilibrium parameters or scientific
configuration interfaces were changed. README additions display four historical
plots with scope qualifications; the existing audit issues remain for Phase 3.

Additions: acceptance criteria, this report, exact transport specification,
direct dependency pins, complete application dependency lock, reproduction
harness, three full-precision profile archives, detailed metrics, input/output
hashes, execution logs and six regenerated figures. Original `results/` is intact.

## Environment and reproducibility commands

Isolated `.venv-phase1`, Python 3.11.1, macOS 26.5.1 arm64; no system site packages.
NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10.3. Transitive application dependencies
are pinned in `requirements-lock.txt`; pip/setuptools versions are in the manifest.
`pip check` executed successfully: "No broken requirements found."
Package installation required network access. pip reported its user cache was
unwritable and disabled caching; this did not prevent installation or checks.

From the repository root, recreate with:

```sh
python3 -m venv .venv-phase1
.venv-phase1/bin/python -m pip install -r docs/phase1/requirements-lock.txt
.venv-phase1/bin/python -m pip check
.venv-phase1/bin/python tools/reproduce_phase1.py --out phase1_runs/NEW_RUN_ID
```

The destination must not already exist. The executed run used
`phase1_runs/2026-11-01-v3-reproduction`. The harness invoked, in sequence:

```sh
.venv-phase1/bin/python hyd.py --out /absolute/path/to/new/run
.venv-phase1/bin/python validation.py --out /absolute/path/to/new/run
.venv-phase1/bin/python compare_h_length.py --out /absolute/path/to/new/run
```

The manifest contains exact argv, working directory, elapsed times, interpreter,
package versions, baseline and execution revisions, settings, timestamps and hashes.
After these scripts, three additional `hyd.run` calls (one per catalyst) generated
full-precision profiles, detailed metrics and sanity checks. These extra calls
are explicitly recorded; they are not solver convergence experiments.
`hyd_contours.py`, `run_all.py` and `--seed-study` were not executed in this phase.
The original comparison script necessarily runs its existing h-only particle
sensitivity; no new transport configurations were implemented.

## Predeclared comparison outcomes

See `acceptance_criteria.md` for tolerances fixed before execution. All four
files match both numerically and **byte-for-byte**, including headers, labels,
flags and line endings. Maximum numeric field difference is zero.

| CSV | Rows | Max absolute difference | Byte equal | Result |
|---|---:|---:|---|---|
| smr_summary.csv | 12 | 0 | Yes | PASS |
| validation_table.csv | 15 | 0 | Yes | PASS |
| validation_equilibrium.csv | 6 | 0 | Yes | PASS |
| h_length_comparison.csv | 6 | 0 | Yes | PASS |

This establishes reproduction at the baseline's saved precision. The original
run's full-precision trajectories and dependency versions are unavailable, so
full-precision trajectory equivalence to the historical run cannot be asserted.
The generated Markdown validation table also matches byte-for-byte.
All six newly generated PNGs differ in bytes from V3; each was visually inspected.
Profiles and comparison trends are consistent with the historical figures;
rendering differs with the recorded environment. This is not pixel equivalence.
The inherited 57-4Q strip-label clipping and the historical annotation/legend
limitations remain; plotting code was deliberately not changed in this phase.

## Physical sanity at baseline conditions

These are executed checks over all saved axial states of PBR, P-R1 and P-R2.
No total-molar-flow conservation condition was imposed.

| Catalyst | Largest relative C/H/O residual | Relative He residual | Largest sum(y)-1 error | Minimum raw flow, mol/s |
|---|---:|---:|---:|---:|
| HMMC | 1.709e-15 | 0 | 2.220e-16 | 0 |
| 57-4Q | 4.556e-15 | 0 | 2.220e-16 | 0 |
| 25-4Q | 2.848e-15 | 0 | 2.220e-16 | 0 |

All states are finite, and all saved raw flows are nonnegative: PASS under the
predeclared thresholds. These checks do not inspect every internal LSODA step.

## Validation metrics reproduced

Experimental errors use five species (CH4,H2,CO,CO2,H2O), full 40 mm P-R2 outlet,
and the digitized experimental-points CSV. MAE and RMSE are mole fractions;
MAPE is the unweighted mean absolute percentage error across the five species.

| Catalyst | MAE | RMSE | MAPE |
|---|---:|---:|---:|
| HMMC | 0.0012172003 | 0.0015795195 | 6.34061% |
| 57-4Q | 0.0039421170 | 0.0050946515 | 16.57526% |
| 25-4Q | 0.0012633181 | 0.0015550441 | 3.85000% |

Figure 4 published simulation agreement is separately evaluated at supplied
PBR digitized points. Species exclude z<=0.05 mm; temperature excludes z<=0.15 mm.
`validation_metrics.json` preserves per-species RMSE, counts and coordinate ranges.
Worst species RMSE is 0.00148444 (HMMC H2). Temperature has two existing domains:

| Catalyst | All available retained temperature points | RMSE K | Existing comparison additionally z<5 mm, RMSE K |
|---|---:|---:|---:|
| HMMC | 89 | 1.235409 | 1.6 |
| 57-4Q | 34 | 2.889192 | 2.9 |
| 25-4Q | 131 | 2.971301 | 3.5 |

Thus the README's <=3 K statement can describe the full retained domain, while
3.5 K belongs to a different restricted domain. The ambiguity is missing domain
qualification, not proof that the two calculations contradict each other.
`species_curve_RMSE` in the comparison table is a mean of five individual RMSEs,
not a pooled RMSE or the largest per-species error.

HMMC H2 reaches the implemented 1% criterion at 1.59 mm. Maximum species
composition difference from the 853 K equilibrium calculation at 2 mm is
0.000479 (saved precision). That equilibrium calculation shares Kp with the
kinetic model, so this is internal consistency rather than independent validation.

## Scientific specification and unresolved issues

See `transport_specification.md` for exact equations, units and paper page/line
references. Main unresolved issues are Re and Nu characteristic lengths,
gas-versus-bed conductivity in Pr, and constant-versus-temperature-dependent ks.
Table 5 internally labels Re and Pr differently from page 7's definitions.
Its SMR Re=1.6 and Pr=.665 give Nu=3.28312 through Eq. 16, rather than 2.91.
No default physics change follows from this observation.

Paper SHA-256: `d43e8fc412676802b7c1f0f8b3562f989e4d7bbdb0708fa0fa55781e5f1efa0e`.
No alternative transport formulation was implemented, evaluated or selected.
No kinetic fitting or equilibrium-constant alteration occurred.

Other unresolved items: 3 mm display/40 mm sampling labels, fixed benchmark
mixing with custom CLI settings, 80 mL/min validation claims, digitization
coverage/uncertainty descriptions, and output-path behavior. They remain
scheduled for Phase 3. P-R2 curve digitization is unavailable; supplied curve
validation is PBR-only. The reference data have figure-reading uncertainty,
not reported physical measurement uncertainty. The actual assignment PDF
and group-specific submission requirements are still unavailable here.

## Scope boundaries and next phase

No numerical reproduction failures occurred. Solver/tolerance convergence,
independent stiff-solver comparison, all-catalyst seed sensitivity, energy
integral consistency and full contour execution remain unverified (Phase 4).
No claim of complete workflow execution or whole-window experimental validation.

Next: Phase 2, only after user approval, implement the named transport sensitivities
while retaining V3 as default; compare identical Figure 4 domains and treat
Table 5/property interpretations separately. Fit quality alone must not select
physical correctness.
