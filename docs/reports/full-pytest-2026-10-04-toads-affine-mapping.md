# Full-suite evaluation — 2026-10-04, source precomputed affine mapping

## Tested implementation

Commit: `7519f2c` (`Match rigid TOADS registration precomputed affine path`).

This change connected the selected rigid TOADS registration path to its
`precompute=true` affine-coordinate behavior. It also corrected the pyramid
Java reference harness, which had called the algebraically equivalent
center-relative `imageToTemplate` expression and therefore did not exercise
the selected class's actual float32 matrix-first coordinate order. The
regenerated Java reference is checked against the Python coarse-to-fine
registration result at every stage and iteration. The Vandy iteration-boundary
test now also checks registration coefficients at a shared Java transform to
within one float32 ULP.

## Full-suite result

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,159 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
857.53 seconds**. Peak sampled RSS was **15.095 GB**; the 24 GB watchdog did not
terminate the run.

The focused registration group passed **32 tests in 11.81 seconds** before the
full run:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_toads_registration_sampling.py tests/test_toads_vandy_iteration_boundary.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_post_registration_branch.py tests/test_toads_alignment_startup.py
```

## Interpretation and remaining work

This is synthetic regression and selected-source slice evidence. It does not
resolve the approximately `1.93e-8` difference in the transform entering the
first base-TOADS registration update, which still changes the strict line
search branch. It does not establish whole-TOADS or MaCRUISE equivalence and
contains no real-MRI validation. The selected TOADS source-line estimate is
unchanged this cycle: existing Java behavior was mapped more faithfully, but
no additional selected Java method lines were newly ported.
