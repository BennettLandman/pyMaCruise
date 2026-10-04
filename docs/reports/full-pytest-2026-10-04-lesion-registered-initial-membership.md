# Full repository pytest — LesionTOADS registered initial membership

- Tested implementation commit: `107cbcf2dafc69a24e92a84ec361a1af41de82e0`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,163 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 872.21s (14:32)**.
- Peak sampled RSS: **13.805 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for the registered LesionTOADS initial
membership composition on synthetic data. It is not whole-workflow Java
equivalence or real-image validation. See the [P39bm report](lesion-registered-initial-membership-2026-10-04.md)
for the focused source-component and composition evidence.
