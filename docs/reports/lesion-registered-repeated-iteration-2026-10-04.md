# Registered LesionTOADS repeated synthetic iteration

- Integration test: `test_registered_continuation_runs_second_synthetic_iteration`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py tests/test_lesion_shapes.py tests/test_lesion13_finite_grow_thin.py`
- Focused result: **17 passed, 0 failed, 0 skipped in 94.74s**.
- Tested implementation commit: `37b857966fda75c04772541d0c8170b1db1cde68`.

Starting from the registered initial membership and first loop, this test
forces the continuing branch through registration/growth/relation refresh and
then resamples the atlas shape maps from the registration result. It invokes
the existing next-iteration thinning, statistics, GM-membership, and
energy-ratio composition and verifies source stage order, an independent
iteration-two thinning state, native fixture-grid dimensions, and finite
13-class memberships/energy ratio. The companion selected 13-class
Java-harness test and the shape-registration tests provide method-level
comparisons. It then runs second atlas registration and an iteration-limit
terminal classification branch; the connected two-iteration path itself is
a synthetic Python integration.

For this compact test, the fixture's synthetic shape arrays are treated as
template-space atlas priors, resampled using the current registration
transform, and the initial constructor distance maps are carried forward.
Multi-iteration convergence and later distance-map refresh behavior remain
open. This is not a two-iteration composite Java reference, complete
LesionTOADS equivalence, or real-image validation. Selected base-TOADS source
line delta: **0**.

Full repository regression evidence is in the [P39bq report](full-pytest-2026-10-04-lesion-registered-repeated-iteration.md).
