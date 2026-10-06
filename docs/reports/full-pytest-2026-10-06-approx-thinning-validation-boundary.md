# Full pytest — approximate queue validation boundary — 2026-10-06

## Tested commit

`9f44477` — `Validate approximate thinning arrays once per pass`

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q --tb=line
```

## Result

- **1,207 passed, 0 failed, 0 skipped**
- **11 expected RAS-orientation warnings**
- **1,217.29 seconds** (20 minutes, 17.29 seconds)
- **14.615 GB peak sampled RSS**, below the 24 GB watchdog and 40 GB machine limit
- Watchdog did not terminate the run

The suite includes the source-derived 343-event approximate queue trace and a check that full-volume membership/shape validation runs once per pass, not once per event. The validation-boundary optimization adds no selected Java source lines and preserves event outputs. This remains synthetic queue-stage evidence, not complete TOADS or whole-workflow equivalence.
