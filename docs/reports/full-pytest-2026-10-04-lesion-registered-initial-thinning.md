# Full pytest — registered LesionTOADS initial thinning

- Tested implementation commit: `abe5bc8d3afdd8e6104cdca4b947b5ba2c71d11b`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,160 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **887.61 seconds**
- Peak sampled RSS: **13.657 GB**
- Watchdog result: normal exit; the 24 GB limit was not reached

The run uses synthetic fixtures and bundled atlas resources. It does not run
the supplied human MRI. The subsequent focused rerun added explicit dtype
assertions to the new integration test; those assertions passed (21 focused
tests), without changing implementation code.
