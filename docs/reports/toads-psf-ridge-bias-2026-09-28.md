# TOADS PSF ridge-bias comparison

**Date:** 2026-09-28  
**Roadmap:** P18af complete  
**Validated implementation/fixture commit:** `9afa5bb1c24c81475f6cbb7b4f45adc691217eee`

## Selected source behavior

`src/pymacruise/toads_math.py::toads_psf_ridge_bias` reproduces the PSF-bias
precomputation in VandyCRUISE `MultipleToadSegmentation.computeExactMemberships`.
For each interior voxel, it checks whether any of the 27 segmentation labels
differs from the center. At such a boundary, and independently for each
channel, it initializes the bias with the center sample multiplied by
`-(27.0f/26.0f)`, adds all 27 neighborhood samples (including the center) at
`1.0f/26.0f` in source nested-loop order, then multiplies by `psfRidge`. When
enabled, the multiplicative field is applied to each image sample before the
stencil. Non-boundary and outer-grid values remain zero. A zero or negative
ridge yields an all-zero field, matching the source `psfRidge > 0` guard.

The selected `MultipleToadSegmentation.java` source SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
`ToadsPsfBiasReferenceHarness.java` transcribes the selected `isBoundary` and
PSF formula excerpts into a deterministic Java harness.

## Synthetic fixture and results

The committed fixture in `tests/data/synthetic/toads_psf_bias/` contains a
9×8×7×2 float32 image, same-grid multiplicative field, and 9×8×7 integer
segmentation as 1 mm isotropic RAS `.nii.gz`. It includes a planar interface,
a small island, uniform regions, and the image border. All three arrays are
checked against the recipe in the committed generator.

Every output value matches the Java harness exactly after float32 conversion
for both field-disabled and field-enabled cases, including the zero-filled
outer faces. Source float32 accumulation can leave a small residual for a
constant-valued boundary neighborhood; the port intentionally preserves that
behavior instead of replacing the stencil with an algebraically simplified
zero. Seven focused tests pass, including zero ridge, non-boundary, border,
and invalid-input cases.

## Scope limits

This validates only the PSF-bias primitive and its boundary decision. It does
not validate the full TOADS membership update, smoothing, shape priors,
normalization, topology-gated label evolution, atlas registration, or final
segmentation equivalence. The committed reference is a source-formula excerpt,
not a full Java segmentation output.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **278 passed in 53.28 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `9afa5bb1c24c81475f6cbb7b4f45adc691217eee`
- Run date: 2026-09-28
