# Full pytest evidence — 2026-09-29 LesionTOADS registration voxel

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **401 passed in 56.22 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `74e15e7a094a71dc14c1784d5f3bf0e1ec748d31` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39h's committed 6³ synthetic case composes LesionTOADS class-wise prior
sampling, template-space derivatives, the rigid coordinate Jacobian, and the
shared registration coefficient accumulator. Priors, gradients, Jacobian,
cost/norm, parameter gradient, and diagonal Hessian match the selected-Java
source-formula excerpt exactly. Class 0 is skipped despite a nonzero map;
zero-prior and class-bound-excluded terms contribute nothing. The fixture and
regeneration instructions are in
`tests/data/synthetic/lesion_registration_voxel/README.md`.

The harness follows selected method-body formulas and does not instantiate
`LesionToadAtlas`. This checks one voxel only. Volume scan order, coefficient
normalization, line search/transform optimization, optimized bounds,
multiple/deformable transforms, and lesion segmentation remain open. The full
run used only synthetic cases; the two warnings are expected orientation
checks from Dura topology tests.
