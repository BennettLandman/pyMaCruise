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
