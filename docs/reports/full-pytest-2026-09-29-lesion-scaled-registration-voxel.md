# Full pytest evidence — 2026-09-29 LesionTOADS scaled registration voxel

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **404 passed in 52.36 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `38db63cdc9ad85e0cb1b3097c4aec6c0af244842` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39j composes one LesionTOADS registration voxel at pyramid scales 2, 4, and
8. The selected-Java formula harness checks mapped template coordinates,
interpolated priors and spatial derivatives, the scaled rigid coordinate
Jacobian, and the voxel's cost, norm, parameter gradient, and diagonal Hessian.
Cases cover multiple registered classes, class-0 exclusion despite a nonzero
map, zero priors, and fractional points outside the template. The formula
fixture and regeneration instructions are in
`tests/data/synthetic/lesion_scaled_registration_voxel/README.md`.

This is a formula-level single-voxel comparison, not a full LesionToadAtlas
run. It does not validate volume traversal/reduction order, line search,
registration optimization, optimized bounds, or lesion segmentation. The full
suite used synthetic cases; the two RAS warnings are expected orientation
checks from Dura topology tests.
