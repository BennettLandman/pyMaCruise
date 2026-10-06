# Full pytest — capped Vandy TOADS output comparison — 2026-10-06

## Tested commit

`84fbba9` — `Validate capped Vandy TOADS output against Java`

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q --tb=line
```

## Result

- **1,203 passed, 0 failed, 0 skipped**
- **11 expected RAS-orientation warnings**
- **1,221.75 seconds** (20 minutes, 21.75 seconds)
- **13.721 GB peak sampled RSS**, below the 24 GB watchdog
- Watchdog did not terminate the run

This regression run includes the 12-iteration synthetic VandyCRUISE output comparison. It is not a whole-workflow source-equivalence or real-MRI acceptance result. Selected TOADS source-line representation is unchanged at an estimated **3,467–4,867 of 13,651 lines (25–36%)**; the change adds no newly mapped Java source methods.
