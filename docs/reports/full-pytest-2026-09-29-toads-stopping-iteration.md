# P18ce full-suite evidence: selected TOADS stopping transition

- **Implementation commit:** `6a91c8d30197c1e3d57409ecd5d0a4c97758b748`
- **Command:** `.venv/bin/pytest -q`
- **Result:** 588 passed, 2 expected RAS reorientation warnings, in 77.38 seconds.
- **Focused check:** `.venv/bin/pytest -q tests/test_toads_stopping_iteration.py` — 3 passed in 0.49 seconds.
- **Composition:** `toads_stopping_iteration` connects selected post-membership control with optional final queue-to-exhaustion growing and final template-label classification. The test asserts source hashes for the selected TOADS segmentation and AlgorithmMultipleToads sources. With `doGrowing=true`, the Java-referenced folded 9³ queue-run fixture matches the composed segmentation and skeleton; with `doGrowing=false`, growth is skipped. Both transfer class indices using the selected mapping and report changed voxels.
- **Limit:** continuing iterations are rejected by this API because their post-growth relation-state refresh is a separate branch. The composition is not a complete TOADS iteration, and no single Java harness executes the combined wrapper.

The Java-matched growing case is documented in [`tests/data/synthetic/toads_growing_pass/README.md`](../../tests/data/synthetic/toads_growing_pass/README.md); loop-decision references are in [`tests/data/synthetic/toads_iteration_control/README.md`](../../tests/data/synthetic/toads_iteration_control/README.md).
