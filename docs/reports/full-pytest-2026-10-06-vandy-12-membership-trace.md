# Full pytest — 12-iteration Vandy membership trace — 2026-10-06

## Tested commit

`eb54ed9` — `Pin per-iteration Java TOADS memberships`

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q --tb=line
```

## Result

- **1,203 passed, 0 failed, 0 skipped**
- **11 expected RAS-orientation warnings**
- **1,216.19 seconds** (20 minutes, 16.19 seconds)
- **13.919 GB peak sampled RSS**, below the 24 GB watchdog
- Watchdog did not terminate the run

This run includes the committed 12-iteration selected-Java membership trace and final-classification hash comparison. It remains compact synthetic evidence, not a whole-workflow source-equivalence or real-MRI acceptance result. Selected TOADS source-line representation remains an estimate of **3,467–4,867 of 13,651 lines (25–36%)**; the slice adds no newly mapped Java methods.
