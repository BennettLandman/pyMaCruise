# Full test suite — corrected LesionTOADS first-thinning handoff — 2026-10-04

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested implementation commit: **`45e917de53ad2b5c81be250627576133ea5eff4a`**.

Result: **1,160 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
863.49s**. Peak sampled RSS: **14.847 GB**. The watchdog did not terminate the
run.

This is regression evidence on the corrected synthetic composition where the
post-constructor membership result is distinct from registration memberships
and the recorded wrapper sequence matches the inspected Java source order. It
does not establish a composite Java output match, full LesionTOADS equivalence,
plus-workflow acceptance, or real-MRI validation. The focused composition
evidence and scope are in the [P39zzzzzzbj report](lesion-registered-initial-thinning-2026-10-04.md).
