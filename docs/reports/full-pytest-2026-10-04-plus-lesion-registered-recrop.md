# Full-suite validation: plus registered LesionTOADS recrop

- **Implementation commit:** `b7ed23eff1ee9dad96e5d6c76c0ba2ef029de18d`
- **Command:** `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,177 passed, 0 failed, 0 skipped; 11 expected RAS reorientation warnings.
- **Duration:** 864.37 seconds.
- **Peak sampled RSS:** 14.806 GB; the 24 GB watchdog did not terminate the run.

P41h adds a restartable plus checkpoint that consumes the normalized-image,
crop, and initial-atlas-alignment manifests. It applies the independent
preprocessing transform for the uncrop and second crop, then the optimized
atlas transform for transformed topology and 13 class shape maps. The
committed synthetic atlas is an 11³ patch from the bundled atlas, keeping this
integration test bounded in runtime and memory.

The focused P41e→P41h synthetic/CLI chain passed 9 tests in 14.92 seconds at
1.299 GB peak sampled RSS. P41h outputs compare bit-for-bit with the existing
`run_registered_lesion_recrop` array adapter; resume skips work when outputs
are intact and repairs a corrupted shape-map file.

This is synthetic source-ordered integration evidence, not an iterative
LesionTOADS execution, a composite Java comparison, or validation against a
human MRI. It does not alter the separate base TOADS convergence findings.
