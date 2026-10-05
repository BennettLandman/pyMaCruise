# Full-suite validation: restartable plus registered LesionTOADS branch

- **Implementation commit:** `f95c4aa4de8c8a66f08f5f685b1772be6b006bd3`
- **Command:** `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,184 passed, 0 failed, 0 skipped; 11 expected RAS reorientation warnings.
- **Duration:** 914.26 seconds.
- **Peak sampled RSS:** 15.274 GB; the 24 GB watchdog did not terminate the run.

P41k adds the `plus-lesion-registered-branch` API/CLI checkpoint. From the
registered recrop it executes initial membership and thinning, the next
statistics and membership pass, rigid shape registration, and one continuing
branch. It saves updated memberships, segmentation, skeleton, availability,
ordering, relation map, summary, and a hash-checked RAS manifest. Focused tests
passed 4 tests in 20.78 seconds; persisted arrays compare exactly with the
in-memory branch result. Resume skips valid outputs, settings mismatches are
rejected, and corrupt output is recomputed. The source-faithful zero-variance
case records its nonfinite energy and zero-deviation diagnostics; the explicit
finite mode completes the synthetic branch.

This is synthetic selected-branch integration. It is not a full plus workflow,
composite Java equivalence, or real-image validation. The checkpoint starts at
the registered recrop and recomputes its own first membership/thinning state;
P41i remains independently restartable at that earlier boundary.
