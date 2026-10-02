# P18cyac: TOADS atlas arrays in RAS order

## Scope

Implementation commit: `a1fa361` (`Load TOADS atlas arrays in RAS order`).

This step adds an algorithm-facing atlas loader for RAS voxel-axis arrays while
retaining `load_toads_atlas_data` for historical Java source-index comparisons.
The RAS API records each asset’s original axis codes, warns when it reorients a
non-RAS map, checks finite 3D voxel data, validates a common grid and integer
topology, preserves template spacing from the selected atlas descriptor, and
stacks shapes in descriptor structure order.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_atlas_assets.py tests/test_toads_alignment_startup.py tests/test_toads_registration_sampling.py
```

Result on the exact implementation commit: **32 passed, 0 failed, 0 skipped in
1.49 seconds**. `git diff --check` passed before commit. The RAS-loading test
uses synthetic 2×3×4 LPS NIfTI assets, verifies three user-visible warnings,
exact voxel reversals, source orientation capture, descriptor spacing, topology
values, and shape-map class order.

## Limits

No human imaging data was used. This validates the array loader and its focused
interfaces, not the complete selected Java atlas setup, optimizer output, TOADS
iteration, or workflow-level numerical equivalence. Descriptor prior assembly
and a single callable startup-to-registration adapter remain open integration
work.
