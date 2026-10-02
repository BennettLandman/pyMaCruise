# P18cyai: atlas-backed registration and template refresh

Implementation commit: `4adb642a69011e15a8b91476ae966ac68e3ade9e`.

`register_toads_from_atlas_assets` composes: read the source atlas descriptor;
load RAS topology and class maps; resolve explicit T1_MPRAGE priors; select the
unique `mask` topology structure as background; prepare the masked one-channel
T1 and initial memberships; run the rigid atlas pyramid; restore the normalized
subject image to its native voxel grid; recompute crop bounds with the fitted
transform; recrop; and sample topology onto the entire recropped grid,
including padding. Image arrays remain in RAS voxel-axis order. The selected
background ratio is caller-visible with default zero. Persistence is outside
this adapter result.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_registration_sampling.py
```

Result on the exact implementation commit: **28 passed, 0 failed, 0 skipped in
1.87 seconds**. The compact synthetic descriptor uses two classes and a folded
RAS subject with anisotropic voxel spacing. Assertions check loaded class
order, registration transform, refresh call order, finite recropped images,
and transformed topology shape equality with the recropped image grid. The
refresh component's selected-Java fixture is included in the focused command.

## Limits

The Java-backed component tests do not make this full descriptor-to-refresh
composition an end-to-end Java comparison. No real subject MRI was used. The
composition stops before segmentation construction, post-refresh initial
memberships, single-scale registration, later TOADS iterations, and outputs.
