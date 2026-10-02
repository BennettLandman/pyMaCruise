# TOADS shape-bound scan vectorization

## Change

`compute_toad_shape_bounds` previously visited every input voxel in Python and
called scalar closest-neighbor interpolation for each one. The scan now batches
up to 65,536 voxels at a time with NumPy. It keeps the selected Java traversal
order, float32 affine arithmetic order, interpolation term order, clipping rule,
positive-support test, and inclusive maximum support index. The bounded batches
avoid allocating full-grid coordinate arrays on human-sized volumes.

## Validation

- `./.venv/bin/python -m pytest tests/test_toads_shape_bounds.py tests/test_toads_shape_volume.py tests/test_toads_default_registration_handoff.py tests/test_toads_lesion_atlas_registration_handoff.py tests/test_lesion_registration_handoff.py -q` — **19 passed, 0 failed, 0 skipped in 3.15s**. This includes selected Java bounds fixtures and the larger folded-prior fixture.
- The vector interpolation kernel was compared bit for bit with the existing scalar kernel at 257 deterministic random in/out-of-bounds coordinates.
- `./.venv/bin/python -m pytest tests/test_toads_initial_registration_control.py tests/test_toads_registered_first_iteration_cycle.py tests/test_toads_default_single_scale_registration.py tests/test_lesion13_finite_grow_thin.py -q` — **5 passed, 0 failed, 0 skipped in 147.47s**.
- Re-ran the packaged fixture's 50-iteration `toads-loop` through the source checkout: **91.63s**, `iterations=50`, `changed_classification_voxels=189`. The output is uint8, finite, 13×13×13, RAS, with labels 1–11. It is bit-identical at every voxel to the pre-change packaged-wheel output (`default-loop.nii.gz`; **0 differing voxels**).
- The source-checkout 2-iteration run took **5.73s** and was bit-identical to the pre-change profiled-wheel output (**0 differing voxels**).

These checks establish output parity for the available synthetic fixture and
selected source-reference cases. The repository-wide full suite has not yet
been rerun on this optimization, so the plan item remains Started until the
next full-suite checkpoint.
