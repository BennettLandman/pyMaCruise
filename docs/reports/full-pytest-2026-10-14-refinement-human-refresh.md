# Full pytest — refinement source correction and human BASE refresh

- Tested commit: `f03cdfe2c1299bd929abcd4752a747968be40700`
- Result: **1,577 passed, 4 skipped, 0 failed**
- Duration: **732.76 seconds**
- Expected orientation warnings: **18**
- Peak sampled RSS: **10.197 GB** under the **16 GB** watchdog
- Scope: full project pytest suite, including the corrected
  `segmentation-refine` CLI expectation based on direct MATLAB `tabulate`
  probes.

Command from the repository root:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 45 -- \
  .venv/bin/pytest -q --tb=short
```

An earlier attempt without `PYTHONPATH=.` stopped at collection because
`tests/test_prepare_plus_tgdm_java_inputs.py` imports the repository's
`scripts` namespace. A subsequent complete run with the repository root on
the module path passed as reported above. An intervening full run exposed a
stale test assertion expecting label 1; that assertion was corrected to 44
after the original MATLAB implementation was directly probed. The focused
refinement/CLI group passed **12 tests** before the final full run.

This suite verifies Python tests, not MaCRUISE whole-workflow equivalence.
Human BASE and held-fixed MATLAB refinement evidence is in
[the human refresh report](in-vivo-base-refinement-refresh-2026-10-14.md).
