# LesionTOADS distance-map refresh after thinning

## Implemented behavior

The selected Java `LesionToadSegmentation.computeMemberships()` refreshes
distance maps from the current hard segmentation when `hasLesion` is true:
ventricle first, brainstem second, then gray matter when `GMflag` is true and
its cap is nonzero. The registered Python transitions now perform that same
call sequence after thinning and before the GM-enabled membership calculation.
The exact-name masks and signed-byte distance kernels are reused from the
already Java-compared Python adapters.

P41k writes the refreshed six-channel maps in its RAS NIfTI checkpoint. P41m
loads those channels on resume, refreshes them after each later thinning pass,
uses the refreshed distance arrays for membership, and persists the new six
channels. A zero distance cap carries the prior map state as the Java wrapper
does; `has_lesion=False` skips the refresh. APIs that omit atlas names retain
their fixed-map behavior for reference callers.

## Validation

- Focused selectors: `tests/test_lesion_distance_refresh.py`,
  `tests/test_lesion_registered_first_iteration.py`,
  `tests/test_lesion_registered_repeated_iteration.py`, and
  `tests/test_plus_lesion_registered_branch.py`: **13 passed** in 88.87s in
  the first focused pass; a final helper plus checkpoint rerun passed **11
  tests** in 73.97s.
- The synthetic P41k integration verifies that at least one distance map
  changes from constructor state after thinning, preserves ventricle →
  brainstem → GM call order, and writes maps equal voxelwise to the selected
  source-matched adapters.
- The synthetic P41m checkpoint verifies stage order before membership and
  exact equality between in-memory refreshed maps and the six persisted RAS
  channels. The existing continuation tests exercise resume and repair.
- Individual ventricle, brainstem, GM-name selection, and capped byte-distance
  behavior remain covered by the adapter tests against their Java references.
  This composes those source comparisons; it is not a whole-wrapper Java
  execution or whole LesionTOADS equivalence claim.
- Full suite on implementation commit
  `483fc25ba12a51316e50ec031672a2e82e7748bc`: **1,191 passed, 0 failed, 0
  skipped, 11 expected RAS warnings in 968.47s**. Peak sampled RSS was **15.270
  GB**, below the 24 GB watchdog. See the [full-suite report](full-pytest-2026-10-05-lesion-distance-refresh.md).

No human MRI was used. Natural convergence, complete LesionTOADS equivalence,
complete MaCRUISEplus execution, and real-data acceptance remain open.
