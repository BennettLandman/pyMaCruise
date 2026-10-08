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

All three meshes are closed, connected genus-zero surfaces under both methods (zero boundary and nonmanifold edges; Euler characteristic 2). Despite effectively coincident vertex positions, the mapped face connectivity remains different. Applying the same Python thickness operator to both mesh sets gives mean thickness 3.5587640 mm for JIST and 3.5582182 mm for Lorensen, a delta of **−0.00054585 mm**.

## Interpretation and limits

The available scikit-image alternative preserves geometry closely but does not remove the local connectivity difference or the resulting thickness shift. This is evidence against switching to Lorensen as a resolution of the P37 substitution risk; it does not prove that the difference is material to downstream segmentation acceptance. Distances sample vertices against the other mesh's triangles and are not a continuous Hausdorff certificate. This is a synthetic, selected-source comparison, not full-workflow or real-MRI validation. P37 remains open until the production substitution is either replaced with source-compatible triangulation or accepted with downstream evidence and an explicit decision.
