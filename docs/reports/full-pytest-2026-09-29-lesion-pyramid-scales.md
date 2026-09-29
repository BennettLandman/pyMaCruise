# Full pytest evidence — 2026-09-29 LesionTOADS pyramid sampling

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **403 passed in 61.34 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `bb35469f4849346bdf3e03fd466e7edcee7e5e9d` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39i compares the LesionTOADS sampler with a Java harness using the selected
`ParametricTransform.imageToTemplate(..., scale)` mapping and
`ImageFunctions` interpolation/derivative formulas at scales 2, 4, and 8.
The synthetic case includes fractional mapped coordinates, anisotropic voxel
resolutions, odd dimensions that truncate on pyramid construction, low/high
constant-zero exterior cases, and explicit full-resolution centers. Values
and in-domain derivatives match exactly. The committed formula fixture and
regeneration instructions are in
`tests/data/synthetic/lesion_scaled_shape_sampling/README.md`.

This tests only scaled prior sampling. It does not implement or validate
LesionTOADS volume-wise coefficient accumulation, the coarse-to-fine
optimizer, optimized shape bounds, multiple/deformable transforms, or lesion
segmentation. The full run used only synthetic cases; the two RAS warnings
are expected orientation checks from Dura topology tests.
