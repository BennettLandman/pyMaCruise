# Full repository test run — 2026-10-18

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
    .venv/bin/pytest -q --tb=short
```

Result: **1,704 passed, 4 skipped, 1 failed, 30 warnings in 960.53 seconds**.
Peak sampled RSS was **11.418 GB**; the 16 GB watchdog did not terminate the
run. The sole failure was
`tests/test_source_reference_plus_driver.py::test_plus_source_run_packet_records_atlasem_as_a_distinct_handoff`, whose
assertion still described `ibsr8obj` as selected after the owner corrected the
target to `ibsr11objV3`. The assertion was updated to verify the new target and
the legacy driver's historical `ibsr8obj` reference. The affected file then
passed **2 tests**. The full suite was not rerun after this test-only assertion
correction.

Separate focused workflow validation passed **106 tests, 0 failures, 4
expected oblique-RAS warnings in 188.45 seconds**. This followed the narrow
implementation fixes from an earlier broad run. These results are distinct:
the targeted passing correction does not convert the full-suite run into a
clean full-suite pass.
