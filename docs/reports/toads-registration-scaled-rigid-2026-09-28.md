# TOADS pyramid-scaled rigid mapping

**Date:** 2026-09-28  
**Roadmap:** P18y complete  
**Implementation commit:** `0854326892f347b95b58d7d582f116674293afa3`

## Source behavior

The selected VandyCRUISE
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/ParametricTransform.java`
defines the rigid pyramid mapping in `imageToTemplate` around line 114. Given
coarse-grid voxel `q` and scale `s`, it first forms full-grid physical
coordinates `(q*s - image_center) * image_resolution`; rotates and translates
that position in millimetres; divides each axis by the template resolution,
adds the full-resolution template center, and divides the result by `s`.
`imageToTemplateDerivatives` around line 198 uses the same full-grid physical
coordinates, then divides rotation derivatives and translation derivatives by
template resolution and by `s`.

The selected source copy has SHA-256
`22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`.
`sample_toad_shape_prior` now accepts an optional positive finite `scale` and
explicit full-resolution image/template centers. Its defaults preserve the
native scale-1 behavior. The already-scaled
`rigid_image_to_template_jacobian` kernel is used unchanged.

## Synthetic comparison

`tests/java_reference/ToadsScaledRigidReferenceHarness.java` is a standalone
formula excerpt of the selected source branches, not execution of the upstream
Java class. Its SHA-256 is
`34a30727b2c61ef7b8e7bcaa8cd633605b6e989e7deeda91827f8e9877fc8eca`. It
generates a committed JSON fixture for scales 1, 2, and 4 using an odd
9×8×7 input grid, odd 11×10×9 template grid, anisotropic resolutions
`(0.8, 1.2, 2.0)` and `(0.9, 1.1, 1.3)` mm, full-grid centers `(4.5, 4.0,
3.5)` and `(5.5, 5.0, 4.5)`, and a nonzero six-parameter rigid transform.
Python mapped coordinates and all 18 Jacobian entries match the harness
float32 values exactly for each scale. The point is chosen on a valid
downsampled z plane even at scale 4, respecting source truncation of partial
blocks. Existing scale-1 shape-sampler fixtures continue to match exactly.
Invalid scales and malformed/nonfinite centers reject.

## Scope boundary

This validates only the scaled coordinate and Jacobian kernels. It does not
integrate a coarse-to-fine optimizer, compare the complete private Java atlas
registration method, establish alignment accuracy, run TOADS segmentation, or
resolve atlas descriptor spacing. Those remain open in P18 and P34.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **263 passed in 57.24 seconds**
- Pytest: 9.1.1; Python 3.13.2
- Tested implementation commit: `0854326892f347b95b58d7d582f116674293afa3`
- Run date: 2026-09-28
