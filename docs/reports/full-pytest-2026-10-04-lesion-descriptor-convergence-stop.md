# Full repository pytest — LesionTOADS convergence-stop branch

- Tested implementation commit: `7e9677a842f306d6e7dd27b0e56eb25d0b4a1ca7`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 874.45s (14:34)**.
- Peak sampled RSS: **13.384 GB**; the 24 GB watchdog did not terminate pytest.

This repository regression run includes synthetic early-convergence and
iteration-limit exit checks for the reusable LesionTOADS descriptor loop. It
does not establish natural convergence behavior on brain-resolution data,
composite Java equivalence, complete LesionTOADS behavior, or real-image
validation. See the [convergence-exit report](lesion-descriptor-convergence-stop-2026-10-04.md).
