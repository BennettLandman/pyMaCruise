# Full pytest: corrected retained-center implementation

- **Implementation commit:** `b785bff` (`Match retained TOADS registration center`)
- **Command:** `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- ./.venv/bin/pytest -q`
- **Result:** 1,503 passed, 0 failed, 0 skipped
- **Expected orientation warnings:** 11, exercising the RAS reorientation warning paths
- **Elapsed:** 1,764.57 s (29 min 24.57 s)
- **Peak sampled RSS:** 13.787 GB under the 16 GB watchdog; the watchdog did not terminate the process

This full-suite run includes the retained-initial-crop-center correction used
by the current brain-resolution TOADS workflow. The full synthetic workflow
integration and resume audit are documented separately in
[the corrected base workflow report](base-workflow-retained-center-2026-10-09.md).

The suite validates the repository's committed synthetic and source-reference
tests. It does not establish complete MATLAB/JIST workflow equivalence or
real-MRI validation.
