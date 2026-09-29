# Full pytest evidence — 2026-09-29 LesionTOADS prior derivatives

| Field | Value |
|---|---|
| Command | `.venv/bin/pytest` from the `pymacruise/` repository root |
| Result | **399 passed in 55.34 seconds** |
| Warnings | Two expected warnings reorienting synthetic LPS Dura inputs to RAS |
| Pytest | 9.1.1 |
| Python | 3.13.2 |
| Tested implementation commit | `683ad5063549de8fcb56f40ec3cdabaab7505fab` |
| Run date | 2026-09-29 |
| Run scope | Complete configured suite; all collected tests passed |

The run includes the P39g LesionTOADS sampler test, which compares its
template-space gradient to the committed selected-Java
`ImageFunctions.linearInterpolation[XYZ]derivative` fixture. The synthetic
case covers fractional and integer coordinates, the exact upper derivative
boundary, low/high out-of-domain zeros, class-bound exclusion, and suppresses
gradient output when the sampled class prior is zero. The derivative helper is
shared with TOADS and has its own hash-checked Java formula excerpt test.

This validates source-slice derivative and sampling behavior only. It does not
validate LesionTOADS atlas registration, registration optimization, or lesion
segmentation. No clinical or user MRI data was used. The two RAS warnings are
expected synthetic orientation cases from Dura topology tests.
