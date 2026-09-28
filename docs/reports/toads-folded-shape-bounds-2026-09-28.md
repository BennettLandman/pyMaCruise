# Folded TOADS shape support bounds

**Date:** 2026-09-28  
**Roadmap:** P18ab complete  
**Validated implementation/fixture commit:** `ecbebfe5ed27d5c22e85271471dc4e2f22da3e84`

## Source behavior

The selected VandyCRUISE `MultipleToadDeformableAtlas.computeTransformedShapeBoundingBox`
scans every input-grid voxel, maps it through the current rigid transform,
samples the class prior with `ImageFunctions.linearClosestInterpolation`, and
includes voxels whose sampled prior is strictly greater than zero. It records
the smallest and largest included integer indices, with the source maximum
stored as the largest index rather than a half-open endpoint.

The added reference harness invokes the selected Java
`ParametricTransform.precomputeImageToTemplateMatrix` and
`ImageFunctions.linearClosestInterpolation` APIs, then applies that bounding
box scan. Source hashes are embedded in `java_reference.json` and checked by
the test. No source checkout or external runtime is needed to run the Python
tests after generation.

## Synthetic validation

The committed RAS `.nii.gz` prior is 39×35×31 voxels at 1 mm isotropic
resolution. It is a deterministic analytic folded ribbon, not a human image.
It is sampled over a 43×39×35 input grid under three transforms:

| Transform | Positive support | Lower bound | Upper bound |
|---|---:|---|---|
| Identity | 13,681 | (3, 0, 0) | (39, 38, 34) |
| Translation (2.25, −1.5, 1.0) mm | 17,769 | (0, 0, 0) | (37, 38, 34) |
| Rotation (0.11, −0.07, 0.05), translation (−1.25, 1.75, −0.5) mm | 19,233 | (2, 0, 0) | (41, 38, 34) |

Python matches every Java support count and each lower/upper coordinate
exactly, and all three transforms produce distinct bounds. This demonstrates
that the fixed-transform scan recomputes bounds when its transform changes.
It does not validate the registration optimizer, the source per-class update
schedule, or TOADS segmentation.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **268 passed in 51.06 seconds**
- Pytest: 9.1.1; Python 3.13.2
- Tested commit: `ecbebfe5ed27d5c22e85271471dc4e2f22da3e84`
- Run date: 2026-09-28
