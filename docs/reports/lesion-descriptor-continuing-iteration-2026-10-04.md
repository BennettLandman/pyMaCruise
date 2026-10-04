# LesionTOADS selected 13-class continuation transition — 2026-10-04

## Scope

This slice adds one callable transition from a continuing post-registration
LesionTOADS branch through the next thinning pass, class/auxiliary statistics
refresh, GM-enabled 13-class membership update, and source float32 energy ratio.
The adapter uses the selected VandyCRUISE behavior already covered by the
adjacent Java-harness tests. It does not add atlas registration, spatial-map or
distance refresh, correction-field fitting, post-membership registration and
growth, loop stopping, file output, or end-to-end LesionTOADS equivalence.

## Validation

Command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_lesion13_finite_grow_thin.py
```

Result: **1 passed, 0 failed, 0 skipped in 87.36 seconds**. The existing test
first checks the 13-class thinning, statistics, and GM-membership state against
pinned Java reference harnesses. The new callable then reruns that same
continuing transition and is compared bit-for-bit against the manual
Java-referenced seam for the thinned labels, class centroids and deviations,
channel weights, all membership and relation maps, objective, and energy ratio.
Its stage trace is checked as thinning → statistics → membership → ratio. The
implementation and focused regression are on local commit
`5f7fc042ba8394435648d9d3b6f6e701d1d5f4df`; no full-suite result is claimed in
the focused evidence. A subsequent full repository run on the same tested
implementation commit passed **1,159 tests, 0 failed, 0 skipped, and 11
expected RAS warnings in 876.20 seconds**, with **14.697 GB peak sampled RSS**
under the 24 GB watchdog; see the [full-suite report](full-pytest-2026-10-04-lesion-descriptor-continuation.md).

## Interpretation

This confirms that the selected component adapters can be assembled into a
deterministic 13-class continuation call on the committed synthetic case. The
comparison is a composition check whose component outputs are individually
covered by the pinned harnesses; it does not independently validate a full
LesionTOADS iteration or a real-image run.
