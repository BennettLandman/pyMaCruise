# Full pytest after unmasked ANTs N4 profile correction — 2026-10-14

The complete suite passed against the tree at commit
`a6a5410eb86e07d240f1126b5e692e0f6a165f4a`, which includes the N4 adapter/CLI
default correction in implementation commit `d4b08f8a023a869fa117a7cc1afdfcd317c8a82c`.

```text
1,577 passed, 4 skipped, 0 failed, 18 expected orientation warnings
675.53 seconds (11 minutes 15 seconds)
Peak sampled RSS: 15.218 GB under a 16 GB watchdog
```

Command, from the repository root:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 MPLCONFIGDIR=build/test-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q
```

This is full regression validation after correcting the no-mask N4 default.
The selected case-100260 human BASE run and its selected Java TOADS
comparisons are reported separately in
[`in-vivo-base-v3-no-mask-n4-source-profile-2026-10-14.md`](in-vivo-base-v3-no-mask-n4-source-profile-2026-10-14.md).
The full suite does not establish whole-workflow Java/MATLAB equivalence or
resolve the failed comparison to the provenance-limited historical output.
The ported TOADS source-method estimate remains 3,511–4,911 of 13,651 lines
(26–36%); no source-method lines were added in this integration cycle.
