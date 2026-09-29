# Full pytest evidence — 2026-09-29

## P18bb: topology-gated TOADS growing event

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **363 passed in 81.61 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `fa1eb5b3cf0b362c11c80361bb5182efea82ede2` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes `tests/test_toads_growing_queue_step.py`, which matches the
single-event Python implementation against the Java harness on topology-accepted,
topology-rejected, and stale-skeleton synthetic states. This is Java harness
reference evidence for the selected source slice, not end-to-end TOADS
equivalence. The generated inputs contain no subject imaging data.

Subsequent plan, report, and README edits do not modify implementation or tests.

## P18bc: initialized TOADS growing pass

| Field | Value |
|---|---|
| Command | `.venv/bin/python -m pytest -q` from the `pymacruise/` repository root |
| Result | **364 passed in 74.98 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `90bb1424ece7378e206748ce7555489474cd34b7` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

This run includes the initialized 448-event growing pass against the selected
Java/LUT/heap harness on a synthetic folded-band phantom. Event outcomes and
all final arrays match the Java harness reference exactly. This is not an
end-to-end TOADS or atlas-registration validation. Subsequent plan, report,
and README edits do not modify implementation or tests.
