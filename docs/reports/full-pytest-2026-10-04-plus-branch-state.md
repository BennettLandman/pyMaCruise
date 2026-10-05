# Full-suite evidence: persisted plus registered-branch state

- **Implementation commit:** `4760612bdb7fb29d9d2cc7c739589ff1f7c27bf0`
- **Command:** `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,184 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 897.13 seconds.
- **Memory:** Peak sampled RSS 14.758 GB; not terminated by the 24 GB watchdog.
- **Focused result:** `PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_lesion_registered_branch.py` — 4 passed in 19.92 seconds.

The focused synthetic checkpoint test compares persisted membership and diagnostic maps against the in-memory API result, checks six persisted distance channels and RAS geometry, verifies energy-ratio/centroid/shape-bound state, resumes without recomputing, and repairs a corrupted output. This is synthetic file-boundary validation. It does not show source-composite equivalence, natural LesionTOADS convergence, or real-MRI validation.
