# Java comparison for scaled TOADS registration kernels

**Date:** 2026-09-28  
**Roadmap:** P18aa complete  
**Validated test/fixture commit:** `3eb7215847917d8d6dda0fa97bbc295b3107f436`

## Source and harness

The comparison combines the selected VandyCRUISE source methods:

- `MultipleToadDeformableAtlas.computeExactScaledRegistrationEnergy` and
  `computeExactScaledRegistrationCoefficients` (`MultipleToadDeformableAtlas.java`,
  SHA-256 `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`,
  around lines 2143 and 2329).
- `ParametricTransform.imageToTemplate` and
  `imageToTemplateDerivatives` (`ParametricTransform.java`, SHA-256
  `22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`,
  around lines 114 and 198).
- `ImageFunctions.subsample`, closest-edge trilinear interpolation, and
  interpolation derivatives (`ImageFunctions.java`, SHA-256
  `be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892`).

`tests/java_reference/ToadsScaledRegistrationKernelsReferenceHarness.java`
is a standalone Java formula excerpt, SHA-256
`16347a6c4f2e1719c3d35368372a32404392e2959f020d1a01959eda2caac3f8`. It
block-averages the committed full-grid synthetic fields in Java float32
offset-loop order, then runs the selected scaled map, interpolation,
derivative, and exact-membership energy/coefficient formulas. It emits the
committed numerical fixture for scales 2 and 4. This is not execution of the
complete private MIPAV/JIST atlas class.

## Synthetic inputs and results

The committed float32 4D RAS NIfTI inputs are in
`tests/data/synthetic/toads_scaled_registration_kernels/`. They have odd
dimensions (memberships 13×12×11×3, atlas shapes 17×16×15×3), anisotropic
spacing, nonzero rigid parameters, and a zero-support tail in the atlas. The
fourth axis is class. Tests confirm RAS axis codes and every input value.

At both scales, Python objective energy, raw cost, and norm match the Java
fixture within `rtol=2e-7`, `atol=3e-8`. Gradient, diagonal Hessian,
norm-gradient, and norm-Hessian arrays match within `rtol=1e-5`, `atol=3e-8`;
the test caps maximum absolute difference across all reported values at
`1.5e-7`. The scalar Python mapper confirms the Java-reported positive-prior
and zero-prior skip counts exactly at both scales. This records a small
float32 coefficient deviation explicitly rather than subsuming it under the
segmentation or surface acceptance limits.

## Scope limits

This validates source-formula kernels on committed synthetic probability
fields. It does not compare the private Java optimizer's line search and
multilevel final transform, test alignment quality, execute TOADS tissue
evolution or dura removal, or resolve the atlas descriptor/NIfTI spacing
discrepancy. Those remain open in P18 and P34.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **267 passed in 55.55 seconds**
- Pytest: 9.1.1; Python 3.13.2
- Tested commit: `3eb7215847917d8d6dda0fa97bbc295b3107f436`
- Run date: 2026-09-28
