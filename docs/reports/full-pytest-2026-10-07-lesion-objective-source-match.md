# Full pytest: LesionTOADS zero-shape objective source match

- Tested commit: `e774747`
- Result: **1,234 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,675.72s** (27m 55s)
- Peak sampled RSS: **10.831 GB**
- RSS watchdog: **16 GB**, not triggered

Command:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  ./.venv/bin/python -m pytest -q
```

This suite includes the selected Java-excerpt comparison for the supported
LesionTOADS membership/objective pass, including zero-shape skips for OBJECT,
LESION, BLACKHOLE, and OUTMASK. It compares all membership arrays, class and
auxiliary energy maps, denominator, running objective, and final objective by
float32 bit pattern. The focused command also passed **10 tests in 6.10s**:

```sh
./.venv/bin/python -m pytest -q \
  tests/test_lesion_objective_membership_pass.py \
  tests/test_lesion_descriptor_membership.py \
  tests/test_lesion_descriptor_iteration_slice.py \
  tests/test_lesion_continuing_iteration.py
```

The evidence is a synthetic selected-method comparison, not full LesionTOADS
workflow equivalence. The source method's atlas generation, statistics
estimation, registration, and later iterations remain separate P39 work.
