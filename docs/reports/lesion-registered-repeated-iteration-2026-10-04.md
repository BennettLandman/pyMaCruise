# Registered LesionTOADS repeated synthetic iteration

- Integration test: `test_registered_continuation_runs_second_synthetic_iteration`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py tests/test_lesion13_finite_grow_thin.py`
- Focused result: **2 passed, 0 failed, 0 skipped in 91.51s**.
- Tested implementation commit: `d962c50be6b0ecdef0f8d79415117854edced255`.

Starting from the registered initial membership and first loop, this test
forces the continuing branch through registration/growth/relation refresh and
then invokes the existing selected next-iteration thinning, statistics,
GM-membership, and energy-ratio composition. It verifies source stage order,
an independent iteration-two thinning state, native fixture-grid dimensions,
and finite 13-class memberships/energy ratio. The companion selected
13-class Java-harness test provides method-level comparison; the connected
two-iteration path itself is a synthetic Python integration.

For this compact test, iteration two receives the fixture's synthetic
transformed shape values and the initial distance maps. Dynamic shape
resampling after registration and repeated-convergence behavior are not
included. This is not a two-iteration composite Java reference, complete
LesionTOADS equivalence, or real-image validation. Selected base-TOADS source
line delta: **0**.

Full repository regression evidence is in the [P39bq report](full-pytest-2026-10-04-lesion-registered-repeated-iteration.md).
