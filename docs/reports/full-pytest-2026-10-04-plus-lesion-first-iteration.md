# Full-suite validation: plus initial LesionTOADS membership and thinning

- **Implementation commit:** `cb558c3eb9ef1a60865d6ea03bafab7f17f6ba16`
- **Command:** `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,180 passed, 0 failed, 0 skipped; 11 expected RAS reorientation warnings.
- **Duration:** 878.87 seconds.
- **Peak sampled RSS:** 15.375 GB; the 24 GB watchdog did not terminate the run.

P41i adds a file checkpoint that consumes the registered recrop stage and
executes the source-ordered initial standard-13 membership, loop-entry, and
first thinning composition. It saves 17 membership channels and the state
needed to continue, with hashes and a restart manifest. The focused P41e→P41i
chain passed 17 tests in 18.04 seconds at 1.423 GB peak sampled RSS. The CLI
help and routing include source loop settings and the finite zero-variance
switch.

The synthetic source-faithful case records a nonfinite objective and stops:
3,968 nonfinite energy values, no nonfinite membership values, 18 zero class
standard-deviation channels, and 2 zero lesion standard-deviation channels.
The separate finite mode uses limiting similarities for exact zero
denominators and proceeds through first thinning. This is an explicit
stabilized deviation from IEEE source behavior, not a numeric tolerance.

This is synthetic first-iteration integration and regression evidence, not a
complete LesionTOADS run, composite Java segmentation match, or real-MRI
validation. The brain-resolution base-run convergence state is unchanged.
