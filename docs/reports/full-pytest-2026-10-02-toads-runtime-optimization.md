# Full synthetic suite: TOADS inverse-prior optimization checkpoint

- **Command:** `./.venv/bin/pytest -q`
- **Tested implementation commit:** `963794e`
- **Result:** 1,119 passed, 0 failed, 0 skipped; 11 expected RAS orientation warnings
- **Elapsed time:** 1,273.87 seconds (21:13.87)

This suite includes the source-order-preserving vectorization of TOADS
inverse-shape-prior accumulation from `6eb502b`, the restartable synthetic
base-run workflow, stage validation/resume behavior, and the previously
validated synthetic source-comparison tests. It was run after the optimization
was added and its focused comparison passed.

The exact source-reference checks and synthetic 10-iteration trace retained
the prior result hash and terminal ratio after vectorization; see the [TOADS
trace report](synthetic-toads-loop-characterization-2026-10-02.md). The
base-run test still reports both TOADS invocations capped and
`converged=false`. This full suite is synthetic regression evidence, not
whole-workflow Java equivalence or real-MRI validation.

The approximate TOADS method-line estimate at this full-suite checkpoint is
**3,434–4,834 of 13,651 selected Java lines**. The vectorization represents an
already-counted source operation, so it adds no new source-line estimate. This
coarse method-level estimate is not runtime code coverage or workflow
completion.
