# Full pytest evidence — 2026-09-29 LesionTOADS registration volume accumulation

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **406 passed in 60.56 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `0740ef256011c30b225ca1d26ac72251b830d2bc` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

P39k adds a source-ordered reducer for pre-sampled LesionTOADS registration
terms. On deterministic 4×3×5 multi-class fields, raw cost, prior-square norm,
gradient, diagonal Hessian, and normalized measures match the selected-Java
formula harness exactly for a dense scan and an offset-zero step-two scan.
The reverse-order harness control differs at float32 precision, showing that
the source x/y/z/class traversal affects the result. A zero-norm test confirms
the selected Java NaN normalization behavior. Fixture provenance and
regeneration instructions are in
`tests/data/synthetic/lesion_registration_volume/README.md`.

This comparison starts from pre-sampled priors, spatial derivatives, and
coordinate Jacobians. It does not validate atlas interpolation, full pyramid
setup, line search, transform optimization, or lesion segmentation. The full
suite used synthetic cases; the two RAS warnings are expected orientation
checks from Dura topology tests.
