# Full repository pytest — LesionTOADS descriptor loop driver

- Tested implementation commit: `dbb6a64741a5befc67487a62d8ff0c0d46f723e0`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 900.55s (15:00)**.
- Peak sampled RSS: **14.083 GB**; the 24 GB watchdog did not terminate pytest.

This repository regression run includes the reusable registered LesionTOADS
loop driver and its compact synthetic comparison. It does not establish
composite Java equivalence, complete LesionTOADS behavior, or real-image
validation. See the [loop-driver report](lesion-descriptor-loop-driver-2026-10-04.md).
