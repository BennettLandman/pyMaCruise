# Full PyMaCRUISE regression suite — compiled TOADS thinning (2026-10-06)

The tested source tree was committed as
`5e0d23eedebfda1104b7c927638f3395280f893b`.

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.2 --report-interval-seconds 120 -- \
  ./.venv/bin/python -m pytest -q --tb=line
```

Result: **1,202 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,182.59s**. Peak sampled RSS was **14.016 GB**. The watchdog exited normally
without terminating the suite. The full suite was run on the exact source tree
later committed under the hash above; no source changes followed the run before
that commit.

Focused checks before the suite:

```bash
./.venv/bin/pytest -q \
  tests/test_toads_thinning_queue_run.py \
  tests/test_toads_initialized_thinning_pass.py \
  tests/test_toads_vandy_iteration_boundary.py \
  tests/test_toads_registered_first_iteration_cycle.py
```

Result: **6 passed, 0 failed, 0 skipped in 13.40s**.

```bash
./.venv/bin/pytest -q \
  tests/test_toads_default_continue_thinning.py \
  tests/test_toads_multichannel_carried_degree_iteration.py \
  tests/test_toads_initial_statistics_thinning.py \
  tests/test_toads_registered_field_refit.py \
  tests/test_toads_p18czm_nonzero_shape_thinning.py \
  tests/test_toads_multichannel_shape_next_loop_refresh.py \
  tests/test_toads_post_registration_next_thinning.py
```

Result: **15 passed, 0 failed, 0 skipped in 17.37s**. These checks validate
the compiled full pass against the source harness and exercise adjacent queue,
registration, and handoff behavior. The source-faithful Python event-trace path
remains available for diagnostics.

This is synthetic regression and Java-harness evidence. It does not establish
whole-workflow Java/MATLAB equivalence or human-MRI acceptance.
