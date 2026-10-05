# Full pytest — 2026-10-05 plus terminal LesionTOADS checkpoint

- Tested implementation commit: `0ace0ffc330bf733de0a4ebebfb0d9cccb2b61eb`
- Command:
  `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,192 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Pytest time: **1,014.58s** (16m 54s)
- Peak sampled RSS: **14.168 GB**; the 24 GB watchdog did not terminate the run.

The warnings are expected input-reorientation warnings exercised by RAS tests.
This suite validates synthetic and unit behavior; it does not establish whole
Java workflow equivalence or real-MRI acceptance.
