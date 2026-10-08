# Full repository test suite: bundled-atlas startup validation

Date: 2026-10-08  
Tested tree: `5a0adcf`

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q
```

## Result

- **1,502 passed, 0 failed, 0 skipped**
- 11 expected RAS-orientation warnings
- Duration: **1,873.83 s** (31:13)
- Peak sampled RSS: **12.442 GB**; watchdog limit: 16 GB
- Watchdog did not terminate the run

This suite includes the committed larger folded-cortex TOADS registration regression. The separate bundled-atlas startup command also completed at a 2.388 GB peak, as documented in [the atlas comparison report](toads-bundled-atlas-registration-java-comparison-2026-10-08.md); that command is not part of the pytest count.

Selected VandyCRUISE source-line coverage estimate remains **3,511–4,911 of 13,651 lines (about 26–36%)**. This cycle adds source-comparison evidence but no newly mapped Java source lines. The estimate is coarse source-method mapping, not runtime Python line coverage or workflow-completion percentage.
