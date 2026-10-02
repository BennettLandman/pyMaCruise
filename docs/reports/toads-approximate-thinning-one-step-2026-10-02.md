# Approximate TOADS thinning: one queue event

## Scope

Added `process_next_approximate_thinning_item` in `src/pymacruise/toads_approximate_thinning_step.py`. It consumes the sparse state produced by `initialize_toads_approximate_thinning`, pops the next source-compatible max-heap item, selects the proposed class from sparse slot 0, applies the shared source-derived relation/topology transition, and updates the removed class volume from its original sparse membership value.

The adapter maps the sparse class/membership pairs into the shared dense local-transition API. It raises the slot-0 score only for target selection; when a removal is eligible, it replaces the resulting volume update with the source's subtraction from the sparse membership associated with the popped class. This lets the shared transition preserve its already-validated local relation fallback and 18/6 topology gate.

## Evidence

- Synthetic input: folded 9 x 9 x 9 three-class label field with two sparse membership slots. The Java harness reproduces the source initialization, max-heap tie path, and one selected `propagateThinning` event. Slot 0 equals the current label; slot 1 can carry a higher membership value, exercising the approximate-path rule that target selection still uses slot 0.
- Compared outputs: popped priority and coordinates/class, class-volume float bits after the sparse decrement, center label/skeleton/availability, and remaining queue size.
- Reference source: `MultipleToadSegmentation.java`, SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- Java harness: `tests/java_reference/ToadsApproximateThinningOneStepReferenceHarness.java`, SHA-256 `44a21b9b7fd5c35e4e2bde88eb45110363d4295f54e4d983325243d414a3aef3`.
- Focused validation: `./.venv/bin/pytest -q tests/test_toads_approximate_thinning_step.py tests/test_toads_approximate_thinning_startup.py` — **3 passed in 0.42 seconds**.

## Limits

This matches one synthetic queue event; it does not establish full approximate TOADS iteration equivalence or imaging-data performance. To exercise removal on this compact phantom, the test supplies zero minimum class volumes. The production Java constructor usually sets each minimum to `max(100, volumeRatio * initialVolume)`, so realistic minimum-volume control still needs an integrated case. After this single event, the adapter does not enqueue newly eligible six-face neighbors and does not run the heap loop to termination. There is no output label change in this fixture because the source's top sparse class equals the current label; it validates skeleton removal and sparse class-volume decrement.
