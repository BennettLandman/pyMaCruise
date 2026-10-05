# Full-suite validation: plus first registered LesionTOADS branch

- **Implementation commit:** `6a0f45a50f072369ad4ec57049ec374a9f7b98b0`
- **Command:** `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,180 passed, 0 failed, 0 skipped; 11 expected RAS reorientation warnings.
- **Duration:** 907.49 seconds.
- **Peak sampled RSS:** 15.775 GB; the 24 GB watchdog did not terminate the run.

P41j carries the explicit finite zero-variance mode from initial membership
through first post-thinning statistics/membership, rigid atlas registration,
and the continuing growth/relation branch on the compact synthetic plus
recrop. The focused chain passed 6 tests in 28.43 seconds. The same case with
no override retains source-faithful NaN objective and energy ratio, confirming
the stabilization is an explicit option rather than silently changing the
historical result.

This is synthetic API-integration evidence. It does not yet provide a
restartable file boundary for the registered branch, validate later plus
iterations, establish composite Java equivalence, or use a real MRI.
