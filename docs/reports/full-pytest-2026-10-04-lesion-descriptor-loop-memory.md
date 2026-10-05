# Full repository pytest — LesionTOADS loop memory bound

- Tested implementation commit: `30d41c8470b80b9f48051d7273e1b758a37fe9fe`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 907.95s (15:07)**.
- Peak sampled RSS: **13.706 GB**; the 24 GB watchdog did not terminate pytest.

This regression run includes the loop-driver compact-history default and the
explicit detailed-history comparison. The RSS measurement is for the complete
repository suite, not an isolated brain-resolution LesionTOADS run. See the
[loop-memory report](lesion-descriptor-loop-memory-2026-10-04.md).
