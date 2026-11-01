# Phase 3 — Code corrections and reproducibility improvements

Status: **PASS for the Phase 3 scope**, 1 November 2026. No Phase 4 execution.
Model-selection decision: `v3_reproduction` is primary/default for reproduction
of published species profiles; the other configurations remain sensitivities.
Neither this selection nor a lower Figure 4 RMSE proves physical superiority.

## File-by-file changes

| File | Change |
|---|---|
| hyd.py | Corrected display/outlet labels; retained PR2_VALIDATION_MM as an alias of PR2_DISPLAY_MM=3; reject custom fixed-benchmark requests before output creation; isolate reference outputs; correct P-R1 region labels; add axes units and margins; attach run metadata. |
| validation.py | Fixed standalone benchmark; signed_error/absolute_error columns and signed equilibrium difference; context-managed CSV reads; exact coverage CSV; full-precision experimental aggregates and per-curve metrics; complete catalyst marker legend and readable figure annotations; metadata. |
| compare_h_length.py | Label the existing sensitivity as Nu-to-h only with tube Re unchanged; name mean_species_RMSE accurately; add units/catalyst identification; isolated outputs/metadata. |
| hyd_contours.py | Remove unsupported 80 mL/min validation wording; distinguish paper-condition stars from supplied benchmark; readable colorbar precision and recovery-contour caption; metadata and safe main entry point. Numerical point/axial calculations unchanged. |
| run_all.py | Main guard; caller-resolved absolute output parent; separate model/fixed_reference/h_length_sensitivity/contours directories; parent/child command and status manifests; no execution on import. |
| run_metadata.py | Shared fresh-output checks, source/data hashes, configuration/solver/dependency/Git provenance, timestamps, output hashes and honest completion/failure states. No reactor calculations. |
| tools/baseline_regression.py | Explicit schema/label mapping to compare unchanged numerical fields to V3; check added absolute_error values. |
| tools/verify_phase3.py | New non-overwriting default regression from a different caller directory; full-array and CSV numeric comparisons; baseline preservation checks. |
| tools/compare_transport.py | Recognize the explicitly documented Phase 3 schema migrations in default regression. Historical Phase 2 files are not edited. |
| tests/test_transport.py | Preserve transport/core regression checks while allowing presentation/schema corrections. |
| tests/test_reproducibility.py | CLI routing, output safety, metadata success/failure, absolute workflow routing, data coverage and plot annotation/layout tests. |
| requirements.txt | Pin the already-used Phase 1 dependency environment; no dependency upgrades or changes. |
| README.md | Clear separation of profile reproduction, experimental comparison and sensitivity; model-selection rationale, metric domains, coverage/uncertainty, current plots, commands and schema/output migration. |

## Physics preservation

No changes to kinetic constants/factors, equilibrium model D, species/energy
balances, property correlations, transport equations/configuration definitions,
solver, seed, grid or reference data. AST checks against Phase 2 confirm numerical
transport, kinetics, equilibrium, balance, integration and contour point/axial
functions are unchanged. Original global physics assignments remain protected
by the transport tests, with the display-constant rename explicitly exempted and
its unchanged value checked. Historical `results/`, Phase 1/2 evidence and
`v3-baseline` remain preserved.

## Executed verification and numeric differences

Criteria were recorded before execution at `ae393e6`. Implementation was committed
at `bed7921`, with label refinement/context-managed readers/README at `b9cff1e`.
The final scientific regression ran at `b9cff1e` (full revision in its manifest).
Primary evidence: `phase3_runs/2026-11-01-corrections-final/`.
An earlier successful regression before layout refinement is separately retained
at `phase3_runs/2026-11-01-corrections/`; it was not overwritten.

Executed from the isolated Phase 1 environment:

```sh
MPLCONFIGDIR=/tmp/hydrogen-phase3-mpl PYTHONDONTWRITEBYTECODE=1 .venv-phase1/bin/python tools/verify_phase3.py --out phase3_runs/2026-11-01-corrections-final
MPLCONFIGDIR=/tmp/hydrogen-phase3-mpl PYTHONDONTWRITEBYTECODE=1 .venv-phase1/bin/python -m unittest discover -s tests -v
```

The regression executes `hyd.py --validate` from the parent directory with a
caller-relative output path; this invokes the standalone reference benchmark
inside `fixed_reference/`. It separately executes `compare_h_length.py`.
It makes three additional default Python `hyd.run` calls for full-array checks.
Exact commands and return codes are in the parent manifest; each actual stage
records its own configuration and output hashes.

| Check | Result | Difference |
|---|---|---|
| Full PBR/P-R1/P-R2 coordinates, flows and temperatures, all three catalysts | PASS | np.array_equal to Phase 1; zero numerical difference |
| Summary CSV numerical fields | PASS | zero vs V3 |
| Experimental validation CSV original numerical fields/flags | PASS | zero vs V3 |
| Equilibrium CSV numerical fields | PASS | zero vs V3 |
| Historical two-case sensitivity CSV numerical fields | PASS | zero vs V3 |
| Original reference data and results preservation | PASS | byte-identical to baseline |
| Source/output manifest hashes | PASS | recorded hashes match files |

CSV **files are not byte-identical**, because section labels and misleading column
names intentionally changed. The regression permits only the documented mappings;
existing numeric differences must equal zero. The new absolute_error field must
match abs(signed_error). The original next-grid-point sample near 3 mm is retained;
no new interpolation alters saved values.

Experimental MAE/MAPE from full-precision benchmark outputs:

| Catalyst | MAE, mole fraction | MAPE |
|---|---:|---:|
| HMMC | 0.0012172003318882427 | 6.340610215997904% |
| 57-4Q | 0.00394211701957371 | 16.575255528787302% |
| 25-4Q | 0.0012633180779300789 | 3.8500005130538875% |

These are five-species experimental outlet metrics, not simulation-curve RMSE.

## Automated tests, visual QA and failure record

Final suite: **14 tests passed** (`test_results.txt`). Coverage includes seven
Phase 2 transport/regression methods and seven Phase 3 reproducibility methods:
custom --validate rejection before output; protected/nonempty/unique output
handling; caller-relative resolution; successful/failed manifest recording;
workflow absolute stage routing; actual CSV coverage; P-R1/field-strip labels;
and operating-window annotations with an explicitly synthetic plotting fixture.
Successful named and legacy CLI calls are executed, not merely parsed.

The first suite attempt had one test error: patching subprocess.run on the shared
module also intercepted Git metadata reads. The production path was not at fault;
the test was narrowed to the workflow module's subprocess reference. Retained
`test_attempt1.txt` records that error. Attempt 2 passed 14 tests but exposed
unclosed CSV reader warnings; readers were changed to context managers, and the
final suite again passed. No warning was reclassified as a failed numerical test.
Upstream Matplotlib/Pyparsing deprecation warnings remain in the logs; dependencies
were not changed merely to suppress them.

All six regenerated reference figures were visually inspected: default profiles,
three field strips, benchmark overlay and h-length comparison. Section labels
are outside data regions, all catalyst markers appear in the overlay legend,
axes have units, and field-strip labels have margins. Exported strip borders
were checked for clipping; all are white at the edges. The final figures replace
active README embeds; historical baseline images are retained as historical files.

Contour annotation coverage uses a synthetic fixture only. The full numerical
contour workflow and run_all pipeline were **not** executed; workflow routing is
tested with mocked children. This is not a claim of Phase 4 completion.

## Interface changes and reproducibility limits

- CLI defaults create timestamped `runs/` destinations; Python benchmark functions
  with no output argument do likewise. Explicit paths must be empty/new and
  outside protected locations. The old overwrite-prone `.`/`results` defaults
  are intentionally replaced; migration is documented in the README.
- `run_all` separates stage outputs; it no longer writes every stage into one
  shared folder. Relative paths are resolved once against the caller.
- `hyd.py --validate` rejects custom T, flow, catalyst selection and nondefault
  transport before writing. Display-only options do not customize the fixed
  reference. Use independent explicit benchmark runs for custom-case comparisons.
- Manifest COMPLETED denotes execution, not physical validation. Failures after
  entering a run context are recorded as FAILED; argument/output refusals before
  starting a run do not create a misleading run directory.
- Low-level plotting/CSV functions do not manage run manifests automatically.
  Their API callers remain responsible for output isolation.
- Exact equality is established in this pinned environment and against available
  Phase 1 trajectories; original V3 full-precision trajectories were unavailable.

## Outstanding issues and next phase

Transport-source/Table 5 ambiguities remain scientific limitations. Visible
figure coverage and digitization uncertainty limit validation; no 80 mL/min
experimental validation or P-R2 curve validation is established. Historical data
and plots retain their old labels by design; corrected active docs supersede them.

Phase 4 still requires full workflow/contour execution, conservation and physical
state checks over the sweep, energy consistency, tolerance/grid sensitivity,
independent stiff solver and all-catalyst H2-seed sensitivity. Also retain the
missing assignment-PDF/report/presentation checklist for Phase 5. No additional
physics or default selection change is proposed for Phase 3 completion.
