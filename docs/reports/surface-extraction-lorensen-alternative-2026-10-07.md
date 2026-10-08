# Lorensen alternative check against selected JIST surface meshes

This synthetic source comparison asks whether scikit-image's classic Lorensen marching-cubes table reduces the triangle-connectivity and downstream-thickness differences seen with the production Lewiner extractor. The candidate and selected Java `IsoSurfaceOnGrid.solveOriginal` consume the same committed 41×41×41 folded RAS level sets. It does not change the production implementation.

## Reproduction

From the `pymacruise` repository root, run:

```bash
./.venv/bin/python scripts/compare_jist_surface_lorensen.py
```

The JSON output, including per-input SHA-256 digests and all measured values, is saved in [`surface-extraction-lorensen-alternative-2026-10-07.json`](surface-extraction-lorensen-alternative-2026-10-07.json). The selected JIST source SHA-256 is `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`.

## Results

| Surface | JIST / Lorensen V/F | Mapped triangle match | Symmetric vertex-to-triangle mean / max (mm) |
|---|---:|---:|---:|
| Inner | 2194 / 4384 | 3018 / 4384 (68.84%) | 1.78e-7 / 8.33e-7 |
| Central | 3062 / 6120 | 4296 / 6120 (70.20%) | 1.73e-7 / 1.61e-6 |
| Outer | 4070 / 8136 | 5624 / 8136 (69.12%) | 1.83e-7 / 1.86e-6 |

All three meshes are closed, connected genus-zero surfaces under both methods (zero boundary and nonmanifold edges; Euler characteristic 2). Under the owner-confirmed P08 metric, the symmetric sampled vertex-to-triangle mean distances (1.73e-7–1.83e-7 mm) pass the 0.01 mm surface-position threshold; the maximum sampled distances (0.83e-6–1.86e-6 mm) are also below the 1 mm Hausdorff threshold. These are the agreed deterministic vertex samples, not a continuous triangle-to-triangle certificate. Despite effectively coincident vertex positions, mapped face connectivity remains different. Applying the same Python thickness operator to both mesh sets gives mean thickness 3.5587640 mm for JIST and 3.5582182 mm for Lorensen, a delta of **−0.00054585 mm**; P08 does not define a thickness acceptance threshold.

## Interpretation and limits

The available scikit-image alternative meets the agreed sampled surface-position tolerances on this phantom, but does not reproduce local face connectivity or eliminate the measured thickness shift. This supports the existing library substitution for the stated geometric tolerance on this case; the connectivity and downstream-measurement differences remain explicit P37 deviations requiring stage-level review. Distances sample vertices against the other mesh's triangles and are not a continuous Hausdorff certificate. This is a synthetic, selected-source comparison, not full-workflow or real-MRI validation. P37 remains open until the production substitution is either replaced with source-compatible triangulation or accepted with downstream evidence and an explicit decision.
