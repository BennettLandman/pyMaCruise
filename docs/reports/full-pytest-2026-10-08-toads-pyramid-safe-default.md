# Full suite: safe-default TOADS pyramid registration

Date: 2026-10-08  
Tested tree: `c62795d` (implementation fix `d42ee16`)

## Command

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q
```

## Result

- **1,501 passed, 0 failed, 0 skipped**
- 11 expected RAS-orientation warnings
- Duration: **1,783.92 s** (29:43)
- Peak sampled RSS: **12.985 GB**; watchdog limit: 16 GB
- Watchdog did not terminate the run

This run validates the complete Python test suite on the committed tree after the coarse-to-fine TOADS registration guard change. It is synthetic and software regression evidence. It does not establish whole-workflow equivalence to MATLAB/JIST or validate real MRI results.
