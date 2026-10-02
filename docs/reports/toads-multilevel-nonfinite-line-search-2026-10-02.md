# TOADS multilevel nonfinite registration trial (2026-10-02)

## Change

`register_toads_atlas_rigid()` now preserves the selected Java line-search
behavior when a trial transform becomes nonfinite. The source computes trial
energy, accepts only when the strict comparison `E > E0` succeeds, rejects a
NaN candidate, and returns a NaN relative change after its bounded retries.
The scale loop then stops because the strict `diff > minEdiff` condition is
false. Python now records `nonfinite_relative_change`, keeps the last accepted
transform, and moves on rather than raising from its coordinate validator.

The source-trial values are pinned in
`tests/data/synthetic/toads_registration_nonfinite_trial/java_reference.json`
to `tests/java_reference/ToadsNonfiniteRegistrationTrialReferenceHarness.java`
by SHA-256. A separate 16³ synthetic case supports all three downsampled
levels (8, 4, 2) and checks finite downstream results after the scale-8 NaN
case. A 13³ CLI test exercises a supported three-level run and verifies finite
uint8 labels on the native RAS grid.

## Evidence

Focused validation command:

```text
./.venv/bin/python -m pytest tests/test_toads_registration_sampling.py tests/test_toads_registration.py tests/test_toads_loop_cli.py::test_toads_loop_cli_runs_supported_multilevel_pyramid_after_nan_trial tests/test_toads_four_level_registration.py tests/test_lesion_output_cli.py tests/test_lesion_output_io.py -q
50 passed, 0 failed, 0 skipped in 42.93s
```

A separate three-level CLI run completed in 43.009s and wrote finite uint8
labels 1–11 on the 13³ RAS input grid. The terminal energy ratio was
0.969482481; the configured threshold was not reached and the iteration cap
was reached.

## Limits

This tests nonfinite line-search handling and supported synthetic pyramid
levels. It does not validate registration convergence, full atlas alignment,
or equivalence on human data. The compact full-loop cap case still reaches its
iteration limit without meeting the energy criterion. The repository-wide
suite has not yet been run on this implementation.
