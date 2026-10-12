# Full pytest suite — LesionTOADS Dura wrapper scale — 2026-10-11

- **Implementation commit:** `3bf4aeb260caaa41cceb3c4668fa7197c3f75fd2`
- **Result:** 1,559 passed, 4 skipped, 0 failed; 11 expected RAS-orientation warnings
- **Duration:** 827.07 seconds (13 minutes 47 seconds)
- **Memory:** 11.211 GB peak sampled RSS under the 16 GB watchdog; watchdog did not terminate the process
- **Command:**

  ```sh
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 16 --interval-seconds 0.2 --report-interval-seconds 60 -- \
    ./.venv/bin/pytest -q
  ```

This suite includes the change that applies the Java
`MedicAlgorithmLesionToads` plugin's ×255 scale to all four Dura wrapper
products, while leaving the normalized low-level generator unchanged. The
synthetic full-plus integration asserts the scaled maps and manifest
provenance. Separate held-fixed P41z downstream evidence shows substantially
improved tissue support and level-set thickness, but ribbon acceptance remains
open. The suite also covers default-on optional regridding and its explicit
opt-out. Passing tests do not establish whole-pipeline source equivalence or
real-human-MRI acceptance.

The selected TOADS method-mapping estimate remains **3,511–4,911 of 13,651
source lines (26–36%)**. This full suite adds no mapped TOADS source methods;
the wrapper correction is in the LesionTOADS plugin boundary. This estimate is
not Python runtime coverage or workflow-completion percentage.
