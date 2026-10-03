# Full pytest evaluation after brain-resolution synthetic integration

## Result

On implementation commit `ad34321c85d4b4e74c0b735a00893ef29073bb7d`, the complete test suite passed:

```text
1138 passed, 11 warnings in 793.64s (0:13:13)
```

There were **0 failed** and **0 skipped** tests. The 11 warnings are expected RAS-reorientation notices exercised by orientation tests. This is synthetic/unit/integration validation; the supplied human MRI was not used.

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.05 --report-interval-seconds 30 -- \
  ./.venv/bin/pytest -q
```

The watchdog sampled a peak of **13.960 GB** before pytest exited successfully, below its 16 GB threshold. The monitor then raced with the already-exited pytest process: a final `proc_pidinfo` call received `ESRCH`, and cleanup attempted to signal the process group. This made the wrapper itself exit with an error even though pytest returned 0 and printed its complete passing summary.

## Watchdog repair

The wrapper now rechecks the child's exit status when RSS lookup reports an OS error, treats a confirmed child exit as normal completion, and avoids signaling an exited process group. A deterministic regression forces the poll/read race. Focused validation on commit `8190a5d0a740abec89ddad160e6d7fd12ad66aae`:

```bash
./.venv/bin/pytest -q tests/test_rss_watchdog.py
```

Result: **2 passed in 0.21s**. The full suite was not rerun after this monitor-only repair; the full-suite pytest result remains tied to `ad34321`.

## Selected source-line accounting

Re-evaluation after this full suite leaves the coarse estimate unchanged: **3,464–4,864 of 13,651 selected TOADS Java lines (about 25–36%)**. This cycle added Numba acceleration for already-represented crop-bound behavior and iteration-level progress reporting, but mapped no new selected Java methods. The estimate counts represented source behavior, not Python statement coverage or end-to-end equivalence; see [the accounting method](toads-source-line-accounting-2026-10-01.md).

The brain-resolution synthetic base route completed all six stages at `(172,216,164)` and resumed with every stage verified/skipped. Its TOADS calls were capped at two iterations. Separate ten-iteration TOADS execution reached terminal energy ratio `0.0216144342` without convergence, peaking at `8.053 GB`. These are synthetic workflow and runtime results, not a source-equivalence or clinical validation claim. Details are in the [brain-resolution runtime report](toads-memory-bounding-2026-10-03.md).
