# Full pytest — approximate TOADS neighbor scheduling — 2026-10-06

## Tested commit

`3f37a7e` — `Port approximate TOADS neighbor scheduling`

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q --tb=line
```

## Result

- **1,204 passed, 0 failed, 0 skipped**
- **11 expected RAS-orientation warnings**
- **1,214.15 seconds** (20 minutes, 14.15 seconds)
- **14.419 GB peak sampled RSS**, below the 24 GB watchdog and 40 GB machine limit
- Watchdog did not terminate the run

The full suite includes the folded three-class approximate thinning event plus face-neighbor scheduling comparison with the selected Java source excerpt. This validates the Python regression suite and one bounded synthetic source comparison; it does not establish complete approximate-path convergence, whole-workflow source equivalence, or real-MRI behavior.
