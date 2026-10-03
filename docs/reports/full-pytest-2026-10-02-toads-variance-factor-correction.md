# Full synthetic suite: TOADS atlas variance-factor correction

- **Command:** `./.venv/bin/pytest -q`
- **Tested implementation commit:** `029255a`
- **Result:** 1,125 passed, 0 failed, 0 skipped; 11 expected RAS-orientation warnings
- **Elapsed time:** 1,265.88 seconds (21:05.88)

The suite includes descriptor parsing for TOADS Intensity Variance Factors,
propagation of those fixed factors through registered membership updates,
the restartable base workflow, and a two-update Java reference comparison.
The comparison checks both exact membership objectives, the source energy
ratio, fixed factor values, and refreshed standard deviations. It validates
the selected update sequence on a small synthetic fixture; it is not a
whole-loop or whole-workflow Java-equivalence result.

The runnable synthetic base command completes all six stages and records
stage output checks, hashes, dependencies, attempt history, and TOADS stop
metadata in `run.json`. A subsequent `--resume` run verifies all six stages
and skips recomputation. Both two-iteration TOADS calls remain capped with
`converged=false`.

The source-line planning estimate is **about 3,464–4,864 of 13,651 selected
Java lines** (about 25–36%). The approximately 30-line increase represents
the selected atlas Intensity Variance Factor accessor and descriptor parser
branches in `MultipleToadDeformableAtlas.java` (lines 270–279 and 1111–1131).
This is a coarse method-level estimate, not executed Java-line coverage or a
workflow-completion percentage.

The full registered TOADS trace remains an open blocker. With the corrected
factor wiring, a 13³ synthetic run capped at 40 iterations took 56.712 seconds
and ended at relative energy ratio 0.0397632, above the configured 1e-4
threshold. Ten-iteration descriptor, uniform-shape, topology-one-hot, and
intensity-scaled probes also did not meet the threshold. These probes are
synthetic diagnostics, not source-equivalence or real-MRI validation. See
[the TOADS variance-factor and convergence report](toads-variance-factor-and-convergence-2026-10-02.md).
