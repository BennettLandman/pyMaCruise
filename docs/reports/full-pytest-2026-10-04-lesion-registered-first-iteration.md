# Full repository pytest — LesionTOADS registered first loop update

- Tested implementation commit: `150e68df478d6248874b596d197f0c2ddfab4c7c`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,164 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 889.67s (14:49)**.
- Peak sampled RSS: **14.024 GB**; the 24 GB watchdog did not terminate pytest.

This records Python repository regression coverage for P39bn's synthetic
composition. It is not a composite Java comparison, complete LesionTOADS
equivalence, or real-image validation. See the [P39bn report](lesion-registered-first-iteration-2026-10-04.md).
