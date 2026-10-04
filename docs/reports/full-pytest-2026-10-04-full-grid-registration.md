# Full-suite evaluation — 2026-10-04, full-grid TOADS registration

## Tested implementation

Commit: `917c15d` (`Match source-wide TOADS registration coefficient scan`).

The post-refresh Python adapter previously restricted exact registration
coefficients to class support bounding boxes. The selected Java method scans
the complete image grid for every registered class. The adapter now uses that
full-grid coefficient path. A new harness calls the selected Java atlas class
and pins its post-refresh coefficients, energy, and output transform against
the same committed synthetic inputs. See the [source-match report](toads-registration-full-grid-source-match-2026-10-04.md).

## Full-suite result

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,159 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
871.73 seconds**. Peak sampled RSS was **14.655 GB**; the 24 GB watchdog did not
terminate the run.

The focused registration group also passed **32 tests in 12.66 seconds**:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_toads_registration_sampling.py tests/test_toads_vandy_iteration_boundary.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_post_registration_branch.py tests/test_toads_alignment_startup.py
```

## Interpretation

This is synthetic source-slice and regression evidence. The post-refresh
transform is now within one ULP of the Java source trace, and the first
registered-loop transform matches bit-for-bit. Post-growth labels, closest
labels, next thinning, and second-iteration memberships still differ. The
result does not establish whole-TOADS or MaCRUISE equivalence and does not
include real-data validation.
