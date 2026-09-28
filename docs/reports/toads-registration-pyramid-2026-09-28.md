# TOADS registration pyramid downsampling

**Date:** 2026-09-28  
**Roadmap:** P18x complete  
**Implementation commit:** `d855faf8e49aab850eb4c47e71ec7fd14702651f`

## Source behavior

The supplied VandyCRUISE source implements 3D and 4D block averaging in
`VandyCRUISE/src/edu/jhmi/rad/medic/libraries/ImageFunctions.java` at line
1946. It creates each output axis with `floor(input_dimension / factor)`, sums
the factor-cubed block in nested x/y/z offset order as float32, divides by a
float32 factor-cubed value, and discards any incomplete boundary blocks.

`registerShapesExactPyramid` in
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`
starts at factor `2**(levels - 1)`, halves to factor 2, and does not run a
full-resolution factor-1 stage. `build_toads_registration_pyramid` reproduces
this schedule and obtains each level directly from the original arrays.

## Numerical validation

`tests/java_reference/ToadsPyramidSubsampleHarness.java` is a standalone
Java transcription of the supplied 4D source method. It generates a
deterministic 9×8×7×3 array and committed JSON output at factors 2 and 4. The
Python `source_block_mean_subsample` result matches every output float32 value
exactly. Factor 4 yields 2×2×1×3, confirming floor truncation on all axes.
Additional tests check coarse-to-fine schedules and dimensions for levels
1–5, array/channel shape agreement, and malformed input rejection. The fixture
is numeric JSON; it is not a medical image.

This step does not port scale-dependent `ParametricTransform.imageToTemplate`
coordinates or run the rigid optimizer at multiple resolutions. Those remain
the next registration integration steps.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **259 passed in 50.60 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested commit: `d855faf8e49aab850eb4c47e71ec7fd14702651f`
- Run date: 2026-09-28
