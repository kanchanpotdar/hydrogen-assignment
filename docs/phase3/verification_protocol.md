# Phase 3 protocol, before execution

No transport equations, kinetic/equilibrium parameters, balances or solver changes.
V3 remains the default; alternatives remain sensitivities. Historical results,
reference data and Phase 1/2 evidence are protected.

Accept default full trajectories against Phase 1: rtol=0 and atol=0 (report exact
array equality); report any deviation as failure. For generated summary CSVs,
allow only explicit section-label changes; all numerical fields must equal V3.
Validation table schema may replace misleading signed-error names with explicit
signed_error and absolute_error. Existing numerical fields must match V3 exactly;
absolute_error must equal abs(signed_error). Equilibrium differences become signed
in name only. The comparison table's species metric becomes mean_species_RMSE;
values must match V3 exactly. Corrected labels/schema do not imply byte identity.

Test invalid custom --validate options fail before writing; named and legacy
interfaces continue to work. Verify resolved caller-relative output directories,
protected-output refusal, unique defaults, child-run paths, completion/failure
manifests and output hashes. Plot tests check annotations and label visibility;
visually review new reference figures. Contour annotation tests may use explicitly
identified synthetic fixtures, not claim full numerical contour execution.
Full contour workflow, convergence, seed and energy checks remain Phase 4.
