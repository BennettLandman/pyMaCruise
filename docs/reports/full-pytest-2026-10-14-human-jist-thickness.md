# Full pytest after human JIST surface and thickness comparisons — 2026-10-14

The complete suite passed against committed tree
`e7ff7dde7fb1b510dd71192d0a0864503f6f3b9e`. This commit contains the human
JIST surface/thickness comparison scripts and P37d/P37e evidence rows.

```text
1,577 passed, 4 skipped, 0 failed, 18 expected orientation warnings
768.89 seconds (12 minutes 48 seconds)
Peak sampled RSS: 11.869 GB under a 16 GB watchdog
```

Command, from the repository root:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 MPLCONFIGDIR=build/test-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q
```

This validates the repository test suite on the recorded commit. The separate
human comparisons are held-fixed Java-source stage evidence; they do not
establish whole-workflow Java/MATLAB equivalence or human anatomical acceptance.
No additional TOADS Java method lines were mapped in this cycle; the selected
TOADS source-method estimate remains **3,511–4,911 of 13,651 lines (26–36%)**.
