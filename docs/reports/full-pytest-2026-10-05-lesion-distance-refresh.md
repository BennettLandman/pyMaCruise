# Full pytest — 2026-10-05 LesionTOADS distance refresh

- Tested implementation commit: `483fc25ba12a51316e50ec031672a2e82e7748bc`
- Command:
  `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,191 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Pytest time: **968.47s** (16m 8s)
- Peak sampled RSS: **15.270 GB**; the 24 GB watchdog did not terminate the run.

The 11 warnings are the intentional input-orientation warnings exercised by the
RAS reorientation tests. This validates the Python test suite on synthetic and
unit-level data; it does not establish complete Java workflow equivalence or
real-MRI acceptance.
