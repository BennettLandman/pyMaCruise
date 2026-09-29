# Full pytest evidence — 2026-09-29 scaled LesionTOADS volume coefficients

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **407 passed in 52.44 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `590bf3c4004fd7013ce6ae540bf8cddd4bf4be3f` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39l composes scaled shape-map sampling, selected rigid coordinate Jacobians,
and P39k source-ordered volume accumulation at fixed transforms. At pyramid
scales 2 and 4, the committed synthetic multi-class cases match the selected
Java `ParametricTransform` and `ImageFunctions` formula harness exactly for
raw cost/norm/gradient/Hessian and their source-normalized outputs. The cases
include odd image dimensions, anisotropic resolutions, class-0 exclusion,
multiple registered classes, zero-default out-of-template samples, and
offset-zero step-two traversal. Fixture provenance and regeneration steps are
in `tests/data/synthetic/lesion_scaled_registration_volume/README.md`.

This is a selected-method, fixed-transform formula comparison, not a complete
`LesionToadAtlas` run. It does not validate pyramid creation, line search,
transform updates, optimized atlas bounds, or segmentation. The full suite
used synthetic cases; the two RAS warnings are expected orientation checks
from Dura topology tests.
