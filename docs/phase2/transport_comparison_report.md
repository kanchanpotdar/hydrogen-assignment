# Phase 2 — Transport sensitivity and scientific comparison

Status: **PASS for Phase 2 scope**. All numerical evidence is in
`phase2_runs/2026-11-01-transport-attempt02/`. No Phase 3 work has begun.

## Implementation and compatibility

The default remains `v3_reproduction`. `Reactor.transport_config` is appended to
the dataclass, preserving existing positional parameters and `h_length` behavior.
A non-None named configuration takes precedence over the legacy Python field;
`effective_transport()` always reports the actual effective choices. CLI defaults
remain equivalent to the original, and conflicting explicit options are rejected.

| Configuration | Re length | Nu-to-h length | Pr conductivity |
|---|---|---|---|
| v3_reproduction | tube | tube | gas |
| v3_h_particle_only | tube | particle | gas |
| particle_lengths_gas_pr | particle | particle | gas |
| printed_equations | particle | particle | bed |

All retain V3's solid-conductivity fit and wall conductivity. `printed_equations`
refers to the dimensionless-group definitions, not every statement in the paper:
it does not implement the separate constant-ks interpretation.

`transport_from_properties()` isolates the closure using SI properties. It returns
U plus Re, Pr, Nu, h, gas/solid/bed conductivities, characteristic lengths, and
configuration. Kinetics, equilibrium constants, governing balances, downstream
heat transfer, seed, numerical tolerances and all original global assignments
are unchanged (verified by AST comparison against the baseline).

Examples, from the repository root:

```sh
.venv-phase1/bin/python hyd.py --transport-config printed_equations --out NEW_OUTPUT
.venv-phase1/bin/python hyd.py --h-length particle --out OTHER_NEW_OUTPUT
.venv-phase1/bin/python tools/compare_transport.py --out phase2_runs/NEW_RUN_ID
MPLCONFIGDIR=/tmp/hydrogen-phase2-test-mpl PYTHONDONTWRITEBYTECODE=1 .venv-phase1/bin/python -m unittest discover -s tests -v
```

The first two numerical interfaces are tested against Python calls. Named
alternatives with `--validate` are rejected because the old validator independently
reruns the default case. Broader custom-case CLI/validation repair remains Phase 3.
No alternative support was added to the old contour runner in this phase; use
the explicit Phase 2 comparison runner for the scientific comparisons.

## Common comparison conditions and evidence

Figure 4: 853 K, 120 mL/min RTP, 5/20/75 CH4/H2O/He, three catalysts, original
solver/seed settings. All configurations use identical digitized point sets;
no interpolation to a common fabricated reference dataset. Per-curve metrics
include sample counts and min/max coordinates. Species use z>0.05 mm;
temperature uses z>0.15 mm, with an additional z<5 mm inlet metric separately.
Experimental MAE/RMSE/MAPE use CH4,H2,CO,CO2,H2O at the full 40 mm outlet.
The mean species RMSE is a mean of five RMSEs, not pooled or worst-species RMSE.

The protocol was committed before evaluation (`330b3ce`). Numerical implementation
was committed at `ab30658`; a comparison-export fix was committed at `86156be`.
The successful comparison ran at `86156bec3cd68b3acb0e1e8b2f432305429410a8`.
Python/dependencies are the isolated Phase 1 environment: Python 3.11.1,
NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10.3. Full dependency pins remain in
`docs/phase1/requirements-lock.txt`. Commands, source/data hashes, configurations,
units/domains, logs, output hashes and full profiles are retained.

## Numerical findings

### Figure 4 agreement and experimental error

| Configuration | Catalyst | Mean species RMSE | T RMSE K, full / inlet | Experiment MAE | Experiment MAPE |
|---|---|---:|---:|---:|---:|
| v3_reproduction | HMMC | 0.0006296 | 1.235 / 1.639 | 0.0012172 | 6.341% |
| v3_reproduction | 57-4Q | 0.0005633 | 2.889 / 2.889 | 0.0039421 | 16.575% |
| v3_reproduction | 25-4Q | 0.0005136 | 2.971 / 3.541 | 0.0012633 | 3.850% |
| v3_h_particle_only | HMMC | 0.0024785 | 12.587 / 16.961 | 0.0012172 | 6.341% |
| v3_h_particle_only | 57-4Q | 0.0015132 | 6.793 / 6.793 | 0.0031413 | 12.572% |
| v3_h_particle_only | 25-4Q | 0.0011583 | 1.263 / 1.701 | 0.0021494 | 4.647% |
| particle_lengths_gas_pr | HMMC | 0.0024374 | 12.401 / 16.711 | 0.0012172 | 6.341% |
| particle_lengths_gas_pr | 57-4Q | 0.0014913 | 6.661 / 6.661 | 0.0031508 | 12.620% |
| particle_lengths_gas_pr | 25-4Q | 0.0011427 | 1.243 / 1.652 | 0.0021385 | 4.637% |
| printed_equations | HMMC | 0.0024347 | 12.390 / 16.695 | 0.0012172 | 6.341% |
| printed_equations | 57-4Q | 0.0014899 | 6.653 / 6.653 | 0.0031514 | 12.623% |
| printed_equations | 25-4Q | 0.0011417 | 1.241 / 1.649 | 0.0021378 | 4.637% |

### Physical predictions and differences from V3

Temperatures are minima on the saved PBR grid, not independently converged minima.

| Configuration | Catalyst | CH4 conversion % | Change, percentage points | Tmin K | Change K | Outlet H2 change, mole fraction |
|---|---|---:|---:|---:|---:|---:|
| v3_reproduction | HMMC | 95.4139 | +0.00000 | 743.5009 | +0.0000 | +0.00000000 |
| v3_reproduction | 57-4Q | 85.3334 | +0.00000 | 831.4244 | +0.0000 | +0.00000000 |
| v3_reproduction | 25-4Q | 67.8195 | +0.00000 | 837.6375 | +0.0000 | +0.00000000 |
| v3_h_particle_only | HMMC | 95.4139 | -0.00000 | 756.2273 | +12.7264 | -0.00000000 |
| v3_h_particle_only | 57-4Q | 86.6140 | +1.28062 | 840.1786 | +8.7542 | +0.00186643 |
| v3_h_particle_only | 25-4Q | 69.3389 | +1.51947 | 843.7754 | +6.1379 | +0.00235514 |
| particle_lengths_gas_pr | HMMC | 95.4139 | -0.00000 | 755.8345 | +12.3337 | -0.00000000 |
| particle_lengths_gas_pr | 57-4Q | 86.5988 | +1.26540 | 840.0377 | +8.6133 | +0.00184426 |
| particle_lengths_gas_pr | 25-4Q | 69.3201 | +1.50068 | 843.6692 | +6.0317 | +0.00232609 |
| printed_equations | HMMC | 95.4139 | -0.00000 | 755.8097 | +12.3089 | -0.00000000 |
| printed_equations | 57-4Q | 86.5978 | +1.26444 | 840.0288 | +8.6044 | +0.00184285 |
| printed_equations | 25-4Q | 69.3189 | +1.49949 | 843.6625 | +6.0250 | +0.00232424 |

### Inlet transport at 853 K

The inlet properties are shared by all catalysts; subsequent local properties evolve.

| Configuration | Re | Pr | Nu | h W/m2/K | U W/m2/K |
|---|---:|---:|---:|---:|---:|
| v3_reproduction | 3.883075 | 0.504199 | 3.717901 | 861.44 | 559.96 |
| v3_h_particle_only | 3.883075 | 0.504199 | 3.717901 | 25843.22 | 1506.72 |
| particle_lengths_gas_pr | 0.129436 | 0.504199 | 2.746372 | 19090.10 | 1476.27 |
| printed_equations | 0.129436 | 0.093805 | 2.701321 | 18776.95 | 1474.37 |

### Table 5, separate from Figure 4

Table 5 lists 773/1173 K and GHSV=2.4e4/h. Correlated calculations infer
135.7168 mL/min RTP from that GHSV and geometry, rather than use Figure 4 flow.
Tabulated calculations use the paper's listed rho/us/mu/cp/kf/ks directly.
No attempt is made to fit rho, feed or velocity to recover its U.

| Property source | T K | Paper U | V3 U | h-only particle U | particle lengths/gas Pr U | printed equations U |
|---|---:|---:|---:|---:|---:|---:|
| tabulated_properties | 773.0 | 584 | 490.51 | 1487.82 | 1466.32 | 1465.12 |
| tabulated_properties | 1173.0 | 1019 | 697.26 | 1533.81 | 1525.50 | 1525.07 |
| correlated_properties | 773.0 | 584 | 548.06 | 1503.79 | 1467.93 | 1465.61 |
| correlated_properties | 1173.0 | 1019 | 656.81 | 1526.91 | 1503.92 | 1502.48 |

All U values above have units W/m2/K. V3 is closest to the tabulated U in these
cases, but none reproduces the complete table. At 773 K, tabulated properties give
Re=1.59775, Pr=.664926 and Eq.16 Nu=3.28242, while the paper lists 2.91.
At 1173 K they give Re=.874151, Pr=.658791, Nu=3.05303 versus listed 2.83.
The correlated densities (.116775/.076954 kg/m3) differ substantially from the
listed .0474/.0206. Correlated velocities (.226396/.343548 m/s) also differ from
.200/.459. Consequently it would be misleading to attribute all Table 5 mismatch
to the Nu characteristic length or to claim that matching U resolves the inputs.

## Which assumptions matter most?

| Controlled change | Observed effect | Scientific interpretation |
|---|---|---|
| Nu-to-h tube -> particle, Re/Pr unchanged | h multiplied by 30; inlet U 559.96 -> 1506.72; Tmin rises 12.73/8.75/6.14 K for HMMC/57-4Q/25-4Q; commercial conversion rises 1.281/1.519 percentage points | Dominant tested assumption. Much stronger heat supply weakens cooling and changes finite-rate conversion. It is the historical hybrid sensitivity, not the complete particle formulation. |
| Re tube -> particle, with particle Nu-to-h length | Re divided by 30; Nu 3.71790 -> 2.74637; U 1506.72 -> 1476.27; commercial conversion falls only .0152/.0188 percentage points | Secondary prediction effect here. Nu's 2.67 term and large wall resistance limit sensitivity to Re. |
| Pr gas -> bed, particle lengths retained | Pr .50420 -> .09381; U 1476.27 -> 1474.37; commercial conversion falls .000965/.001192 percentage points | Small incremental prediction effect at this inlet, despite large dimensionless Pr change. This does not establish that the definitions are interchangeable. |

In the particle-length cases U approaches the wall-resistance ceiling
kw/s=1600 W/m2/K. That explains why further Re/Pr changes have limited leverage
on U here. This ranking is conditional on this geometry, properties and case;
it is not a universal ranking over operating conditions.

HMMC reaches effectively the same equilibrium outlet under all choices (outlet
H2 changes below 6e-13). Its axial profiles and sampled thermal dip still change
substantially, so outlet agreement alone cannot identify the correct transport.

All alternatives worsen published species-profile agreement. They improve the
57-4Q experimental outlet MAE (about .00394 -> .00314/.00315), yet worsen the
25-4Q experimental MAE (about .00126 -> .00214/.00215). Temperature agreement
improves for 25-4Q but worsens for HMMC and 57-4Q. There is no universal winner.

## Physical defensibility and recommendation

- Keep V3 as the reproduction default. Its Table 5-style Re/Pr choices and closer
  Figure 4 agreement support its role as an explicit reproduction interpretation,
  not proof that its heat-transfer closure is physically definitive.
- `particle_lengths_gas_pr` is a more internally consistent particle-length
  alternative than the historical h-only switch; gas Pr retains a fluid-property
  definition. It remains an interpretation requiring source confirmation.
- `printed_equations` follows Costamagna Eqs.11–13 and their particle-length text
  more literally, but conflicts with Table 5's definitions. Its name is not a
  physical-correctness certificate.
- The original Singhal correlation concerns gas-particle heat transfer. Applying
  it through an effective bed conductivity in a wall/furnace resistance model
  requires care. The author institution's abstract supports that context:
  https://www.sintef.no/en/publications/publication/0198cc3fa21e-f7472bdc-19b2-4010-9c8d-c1b48aae70c0/
  The publisher full text and linked repository could not be retrieved (403/cache
  failure). This phase does not claim independent verification of the original
  correlation's exact definitions, calibration range or applicability here.

No alternative is established as physically superior by this evidence. Reconcile
source definitions and intended wall-versus-particle closure before proposing a
new default. Constant ks=1.6 is not investigated: it would add another material
assumption without resolving those ambiguities, and the agreed four cases retain
V3's ks fit. No coefficient tuning or new governing equations are justified here.

## Verification, failures and limitations

Seven unittest methods passed: per-configuration equations/effective parameters;
dimensional scaling and wall resistance; legacy aliases and explicit precedence;
CLI choices/conflicts; successful named/legacy CLI runs; default full profiles
and summary regression; AST preservation of the scientific core and globals.
See `test_results.txt`. The full-profile comparison uses the predeclared Phase 2
tolerances against Phase 1; it does not establish full-precision equivalence to
the unavailable original V3 trajectories.

The three original scripts were executed without custom options into a new
legacy-regression directory. All four numerical CSVs match the baseline tag
byte-for-byte (zero differences). All 12 configuration/catalyst state checks pass;
maximum relative elemental error 4.5565e-15, helium residual zero, no nonfinite
or negative saved flows. All three new scientific comparison figures were
visually inspected. No protected reference-data or original output files changed.

The first attempt, `phase2_runs/2026-11-01-transport/`, failed after legacy
regression because the reporting export duplicated a configuration keyword.
Its FAIL manifest and logs are retained. Commit `86156be` fixes that export;
the new attempt completed successfully. The failed attempt is not relabelled.

No full contours, solver convergence, independent stiff solver, or energy-integral
verification were run in Phase 2. No 80 mL/min experimental validation is claimed.
Digitized PBR curves are not original author numerical data; P-R2 curves are
unavailable. Five measured species and incomplete visible temperature coverage
limit physical validation. Existing baseline annotations and general custom CLI
validation issues remain for Phase 3. The historical Phase 1 harness intentionally
refuses modified model sources; reproduce that historical phase at its recorded
revision. Use the Phase 2 runner to verify the current default against V3.

Next phase, subject to user approval: Phase 3 code/documentation corrections and
broader reproducibility improvements. Do not change governing physics there.
