# TOADS initial exact membership through first thinning (P18cyar)

P18cyar composes the selected source sequence from initialized class state
through one exact membership pass and the first thinning transition. It
computes relation, boundary-distance, and closest-label state; uses those
closest labels in the initial deviations/image-variance/similarity/membership
composition; and supplies the resulting memberships plus relation adjacency
to the selected thinning implementation.

The test uses the folded asymmetric 9×8×7 two-channel, 11-class RAS dispatch
fixture, a synthetic assigned-label state, and explicit transformed shapes.
It checks order, relation handoff, membership handoff, thinning outputs against
direct calls to the same components, and input immutability. The relation,
membership, and thinning components have separate selected-Java reference
tests. This integration test is not a combined Java wrapper comparison and
does not demonstrate topology-changing behavior on a large volume.

- Command: `./.venv/bin/pytest -q tests/test_toads_initial_statistics_thinning.py tests/test_toads_initial_statistics_membership.py tests/test_toads_initial_relation_state.py tests/test_toads_default_atlas_dispatch.py tests/test_toads_default_continue_thinning.py`
- Result: **14 passed, 0 failed, 0 skipped in 2.87 seconds**.
- Implementation commit: `ada020ac149e2889ede1195202101bc702adcc28`.
