# Full pytest evaluation — lesion type normalization — 2026-10-11

## Result

The full test suite passed against implementation commit `0929c81f`:

- **1,543 passed, 4 skipped, 0 failed, 11 expected RAS-orientation warnings**
- **818.41 seconds**
- **10.921 GB peak sampled RSS**, under the 16 GB watchdog limit
- Watchdog did not terminate the run

The working tree contained only documentation edits during this run; no
implementation changes were made until after the suite completed. This
evaluation includes the lesion object-type normalization fix and its
regressions. It predates the subsequently identified PLUS fusion-mask wiring
correction, which requires a focused regression and source comparison.

## Reproduction

From the `pymacruise` repository root:

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 \
  MPLCONFIGDIR=build/p41z-case-normalization/full-suite-mpl \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
    .venv/bin/pytest -q --tb=short
```

The watchdog sampled the process throughout the run; its final output was
`child_exit=0 rss_limit_gb=16 peak_sampled_rss_gb=10.921 terminated_for_rss_limit=false`.

## Ported source coverage estimate

The selected TOADS source-method mapping remains **3,511–4,911 of 13,651
lines (26–36%)**. This full-suite run did not add TOADS source-method coverage;
it validates the implementation regression set and is not a workflow-completion
percentage.
