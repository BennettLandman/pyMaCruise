# Full pytest — TOADS JAMA LU correction — 2026-10-06

- Tested tree: `2b1d0c51888725f3d6a79138ac777bab06bb0b0d` (implementation `52b68f5d1b74d723c2a93e76031b16ac2e25e193` plus status/evidence commits)
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- ./.venv/bin/pytest -q`
- Result: **1,218 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,452.34s** (24m 12s)
- Peak sampled RSS: **10.179 GB**; the 24 GB watchdog did not terminate the suite.

This is repository-wide Python regression evidence after adding source-ordered JAMA LU solving for TOADS polynomial field fits and exact compact Vandy loop-trace comparisons. It does not establish base/plus whole-workflow equivalence or real-MRI acceptance. The source-behavior estimate for the eight selected VandyCRUISE TOADS files remains **3,511–4,911 of 13,651 lines (about 26–36%)**: this change represents numerical dependency behavior from JAMA, which is outside that denominator.
