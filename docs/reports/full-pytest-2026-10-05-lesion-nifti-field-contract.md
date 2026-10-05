# Full PyMaCRUISE regression suite — 2026-10-05 LesionTOADS NIfTI field contract

Tested implementation commit: `24007918da04195c8dae2df06835f36694fbd3eb`.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,202 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,232.66s**. Peak sampled RSS: **12.750 GB**. The watchdog exited normally and
did not terminate the run.

The suite includes optional LesionTOADS correction-field persistence and
automatic restore across registered-loop checkpoint segments; workflow and
single-segment CLI routing; rejection of uncompressed correction-field input;
and regression coverage for a child exiting between watchdog polling and RSS
sampling. The focused registered-loop file passed **12 tests in 127.81s** on
implementation commit `abbfc7830a11959cd0cbe97b3255df62d1f9dc48`. The focused
watchdog file passed **2 tests in 0.52s** on the implementation commit above.

This is synthetic Python regression evidence. It does not establish composite
Java equivalence or human-MRI validation.
