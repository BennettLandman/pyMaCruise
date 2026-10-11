# Full-suite result for corrected in-vivo Dura profile

## Suite run

- **Tested implementation commit:** `3bb0634` (`Resume plus pipeline after repaired inputs`)
- **Command:**

  ```sh
  PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
    ./.venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 24 --interval-seconds 0.2 --report-interval-seconds 120 -- \
    ./.venv/bin/pytest -q
  ```

- **Result:** 1,512 passed, 1 failed, 0 skipped, 11 warnings
- **Duration:** 1,137.75 seconds (18 minutes 57 seconds)
- **Peak sampled RSS:** 12.730 GB; the 24 GB watchdog did not terminate the run.
- **Source-line delta:** 0 newly mapped selected Java method lines. This was integration and restart behavior validation, not new algorithm-port coverage.

The sole failure was `tests/test_plus_file_workflow.py::test_plus_file_workflow_runs_resumes_and_repairs_invalid_tail`. Its final assertion expected changing `crop_edge_width` during resume to raise `ValueError`. The current workflow intentionally accepts the changed run signature and uses stage dependency signatures to invalidate and recompute affected checkpoints.

## Corrected assertion

The test now changes the crop width from 5 to 4 and asserts that no stages are reused, the new run completes, the prior signature is retained for provenance, and the manifest records the new setting. Focused verification after that test-only change:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  ./.venv/bin/pytest -q \
  tests/test_plus_file_workflow.py::test_plus_file_workflow_runs_resumes_and_repairs_invalid_tail
```

**Result:** 1 passed in 41.37 seconds.

The full suite has not been rerun after this test-only correction. The 1,512/1 failing result above remains the most recent full-suite evidence; the focused result verifies the corrected behavior and assertion.
