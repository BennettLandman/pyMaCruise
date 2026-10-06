# Full pytest — approximate TOADS growing — 2026-10-06

- Tested implementation commit: `9c7a800bbb7b33b774a3a75fdf6e8fb74e48a350`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- ./.venv/bin/pytest -q`
- Result: **1,209 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,309.68s** (21m 49s)
- Peak sampled RSS: **13.363 GB**; 24 GB watchdog did not terminate the suite.

This repository-wide regression result covers the approximate sparse-growing
implementation and its 480-event source-derived Java fixture. It does not
validate whole TOADS workflow equivalence, approximate-mode convergence, or
real imaging data.
