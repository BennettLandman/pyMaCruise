# LesionTOADS registered third synthetic iteration

- Integration test: `test_registered_continuation_runs_second_synthetic_iteration`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py tests/test_lesion_shapes.py tests/test_lesion13_finite_grow_thin.py`
- Focused result: **17 passed, 0 failed, 0 skipped in 97.51s**.
- Tested implementation commit: `6def02db2281f804efbcc40d506bf3fd7ebcb7c5`.

The RAS 9³ synthetic fixture now carries the registered LesionTOADS state through
three successive membership updates. It refreshes image-space shape values at
the first and second registration transforms, verifies the second registration
continues, checks finite third-iteration class memberships and energy ratio, and
routes the third registration through the configured iteration-limit terminal
classification branch. The fixture keeps constructor distance maps fixed and
uses source-zero-variance handling disabled for these non-degenerate synthetic
inputs.

This is synthetic Python composition over separately source-checked methods.
It is not a composite Java output comparison, an unconstrained convergence
test, complete LesionTOADS validation, or real-image validation. Full suite:
**1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 903.35s**,
peak sampled RSS **14.892 GB** under the 24 GB watchdog. See the
[full-suite report](full-pytest-2026-10-04-lesion-third-iteration.md).
