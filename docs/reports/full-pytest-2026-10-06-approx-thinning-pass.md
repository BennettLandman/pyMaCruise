# Full pytest — approximate TOADS queue pass — 2026-10-06

## Tested commit

`2e16f4b` — `Compose approximate TOADS thinning queue pass`

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q --tb=line
```

## Result

- **1,206 passed, 0 failed, 0 skipped**
- **11 expected RAS-orientation warnings**
- **1,212.71 seconds** (20 minutes, 12.71 seconds)
- **14.427 GB peak sampled RSS**, below the 24 GB watchdog and 40 GB machine limit
- Watchdog did not terminate the run

This suite includes exact event-trace-hash and terminal-state comparison of one bounded synthetic approximate queue pass with a source-derived Java harness. It does not validate all approximate `propagateThinning` branches, the full TOADS loop, whole-workflow equivalence, or real-MRI behavior.
