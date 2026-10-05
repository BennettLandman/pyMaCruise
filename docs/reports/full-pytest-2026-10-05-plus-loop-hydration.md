# Full-suite evidence: hydrated plus registered-loop continuation

- **Tested commit:** `9a854c9`
- **Full command:** `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,184 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 923.59 seconds.
- **Memory:** Peak sampled RSS 15.276 GB; not terminated by the 24 GB watchdog.
- **Focused synthetic test:** `PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_lesion_registered_branch.py::test_registered_branch_checkpoint_completes_and_resumes` — 1 passed in 16.34 seconds on code/test commit `1426d68`.

The focused test hydrates persisted P41k arrays/scalars and checks them against captured in-memory state. It then runs iteration 2 through the registered-loop API and compares segmentation, class memberships, and registration transform against a direct loop call starting from the in-memory first-iteration branch. This is synthetic API-composition evidence. The next branch is not yet persisted, no natural convergence claim is made, and no human MRI was used.
