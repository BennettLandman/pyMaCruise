# Full PyMaCRUISE regression suite — TOADS ordering carry (2026-10-05)

Tested implementation tree, committed as `c8a33033e823d9a2fbce5170ee7ebd984008e801`.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,202 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,187.29s**. Peak sampled RSS was **13.539 GB**. The watchdog exited normally
without terminating the suite.

The full run was executed immediately before commit. After the run, the only
change before commit corrected a stale scalar in the superseded historical
Python-diagnostic block of the Java-reference JSON; this scalar is not used by
the current assertions. The two tests that parse and validate that reference
were rerun after the correction: **2 passed in 12.15s**. The selected Java
comparison and regenerated Python folded-cortex snapshot are included in the
committed tree.

This is synthetic regression evidence. It is not a real-MRI validation, and
not whole-MaCRUISE equivalence.
