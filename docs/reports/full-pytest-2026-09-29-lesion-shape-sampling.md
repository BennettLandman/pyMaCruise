# Full pytest evidence — 2026-09-29 LesionTOADS shape sampling

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **398 passed in 55.53 seconds** |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `9a71937c62edb34c082f09c883624307d0236e9b` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

The run includes eight synthetic comparisons against a Java formula harness
for the selected LesionTOADS rigid shape-sampling wrapper. It verifies mapped
coordinates and exact float32 values for identity/nonzero rigid transforms,
anisotropic resolutions, per-class image bounds, class-specific default
values, and low/high out-of-template coordinates. It also checks direct
untransformed indexing with a RAS NIfTI array and rejects unsupported
transform modes and invalid maps.

The harness mirrors the selected `LesionToadAtlas.getRegularTransformedShape`
wrapper around selected Java transform/interpolation APIs; it does not compile
or instantiate the full atlas class. This validates one fixed rigid sampling
slice, not atlas registration, optimized bounds, deformable transforms, or
LesionTOADS segmentation.

Two expected warnings come from a synthetic LPS-oriented Dura input and paint
mask, which the suite confirms are reoriented to RAS. No clinical or user MRI
data was used.
