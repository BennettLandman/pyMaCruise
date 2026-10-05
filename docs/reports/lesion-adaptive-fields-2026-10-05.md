# Adaptive LesionTOADS field refitting in the registered loop — 2026-10-05

## Implemented call order

`LesionToadsDescriptorLoopOptions` now optionally accepts a maximum correction
polynomial degree and an initial carried degree. When enabled with attached
channel-last correction fields, each continuing iteration runs the selected
degree gate and exact global Chebyshev refit after thinning and before
`computeCentroids`. The refit consumes the previous membership maps, newly
thinned class labels, and current class centroids in source channel order.
Statistics and the next GM-enabled membership pass then use the updated
fields. The loop carries degree/field state to later passes and returns the
final fields and degree for callers to persist. It copies input fields before
refitting. With the feature unset, the loop call path is unchanged.

The selected plus file workflow continues to disable correction by default.
Exposing these controls in restartable workflow configs/manifests and persisting
the fitted fields at checkpoint boundaries remain separate work.

Implementation commit: `2b90545a52ace6a3d518925a380c9c34ccac88a0`.

## Focused validation

```bash
PYTHONPATH=src ./.venv/bin/pytest -q --tb=short \
  tests/test_lesion_registered_repeated_iteration.py::test_registered_loop_driver_matches_two_manual_updates \
  tests/test_lesion_inhomogeneity_transition.py \
  tests/test_lesion_toads_field_data_terms.py \
  tests/test_lesion_toads_fielded_data_energy_membership.py \
  tests/test_lesion_standard13_voxel.py \
  tests/test_lesion_standard13_volume.py \
  tests/test_lesion_post_thinning_gm_membership.py \
  tests/test_lesion_descriptor_iteration_slice.py \
  tests/test_lesion_toads_data_terms.py \
  tests/test_lesion_active_class_energy.py
```

Result: **25 passed, 0 failed, 0 skipped in 31.43s**. The synthetic loop
test enables degree-one fitting, checks refit stage order before centroid
refresh, verifies finite returned field/degree state, and checks the caller's
initial fields remain unchanged. Existing Java-reference tests in the group
cover selected-source fielded energy terms and adaptive gate/refit components.
This is not a composite Java algorithm comparison. Full suite is pending for
this implementation.
