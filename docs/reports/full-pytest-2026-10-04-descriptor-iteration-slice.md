# Full pytest after LesionTOADS descriptor iteration slice

Command:

```text
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested commit: `3aadd8b18e3d80fdb90a5f5fe2ed9a9faa0bf78d`.

Result: **1,150 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
846.38 seconds (14:06)**. The watchdog exited 0, sampled peak RSS of
**12.444 GB**, and did not terminate the process (24 GB cap).

This suite exercised the newly composed raw-descriptor-energy to
membership-update slice and all repository synthetic regressions. No human
MRI data was run. The new Java harness validates a single interior voxel with
fixed centroids, weights, shapes, similarities, and direct relations. It is
not iterative LesionTOADS equivalence, full workflow equivalence, or
real-data validation. No selected base-TOADS source lines were newly counted;
the current estimate remains 3,464–4,864 of 13,651 (25–36%).
