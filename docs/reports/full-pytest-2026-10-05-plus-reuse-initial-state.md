# Full-suite evidence: reuse persisted P41i initial state

- **Tested commit:** `2b4cf5492ed9117b99c301b4248a3ab1f38b7b5d`
- **Full command:** `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,186 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 966.37 seconds.
- **Memory:** Peak sampled RSS 15.406 GB; the 24 GB watchdog did not terminate the run.
- **Focused chain:** `PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_lesion_first_iteration.py tests/test_plus_lesion_registered_branch.py` — 8 passed in 64.59 seconds.

P41i now records membership energy maps, six initial distance channels, relation IDs/table, and first-thinning class volumes. P41k accepts the P41i manifest, validates that it matches the recrop and selected settings, and reconstructs the in-memory initial state without rerunning initial membership/thinning. The persisted-state branch matches the direct in-memory transition membership, registration transform, and segmentation on the compact synthetic RAS case. The original recrop-start path remains available. This is synthetic selected-branch integration, not full plus or source-composite equivalence, and does not use human MRI data.
