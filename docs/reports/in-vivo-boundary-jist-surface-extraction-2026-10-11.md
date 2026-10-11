# In-vivo selected JIST surface extraction comparison

This compares the saved human base-tail `inner`, `central`, and `outer` level sets from `build/invivo-base-source-v3-human-boundary-2026-10-11` against both the production Python VTK meshes and a direct replay of the selected Java `IsoSurfaceOnGrid.solveOriginal` source on identical float32 values. This corrected tail supplies the cleaned MaSeg boundary to ACE mode 2. It holds the level sets fixed; it is not a full workflow run or evidence of upstream TOADS, ACE, refinement, thickness, or clinical equivalence.

Selected source SHA-256: `a323ac5eed8b66ad62a68cd93a6bc7bcf01281453fc6541ee5c4ccc529d40c38`. Wrapper SHA-256: `3f7b54e099f2550f42a572e96070642b8124a93ab1b45eca1aa544ecec1758b4`.

The Java runner uses the selected source body with small JIST container/progress/normal-generator stubs. The field is cropped around all negative samples with a two-voxel positive halo; each crop boundary is verified strictly positive, so the crop does not truncate a zero crossing. Java heap is limited to 2 GiB. Java vertices are translated by the crop start and transformed with the source NIfTI affine to RAS millimetres. The Python meshes are the saved production RAS-mm VTK outputs from the same level sets.

## Geometry and topology

| Level set | Java V/F/components/boundary/nonmanifold/χ | Python V/F/components/boundary/nonmanifold/χ | Symmetric vertex-to-triangle mean (mm) | Sampled Hausdorff (mm) | User limits |
|---|---:|---:|---:|---:|---|
| inner | 234234/468468/143/0/0/0 | 234234/468468/143/0/0/0 | 0.00000012 | 0.00034527 | pass |
| central | 235178/470356/143/0/0/0 | 235178/470356/143/0/0/0 | 0.00000013 | 0.00034527 | pass |
| outer | 303098/606196/143/0/0/0 | 303098/606196/143/0/0/0 | 0.00000014 | 0.00048828 | pass |

The reported symmetric surface statistic samples every vertex of each mesh and queries its closest point on the opposite triangle mesh using PyMaCRUISE's selected-JIST triangle-tree distance implementation. The maximum is a sampled vertex-to-triangle Hausdorff estimate, not a formal continuous triangle-to-triangle bound. Topology counts include all connected components, making small islands visible instead of silently discarding them.

## Scope and limitations

- This directly replays the Java `IsoSurfaceOnGrid` algorithm source, not the complete JIST runtime, VTK serializer, or historical full MaCRUISE pipeline.
- Both extractors consume fixed saved human fields. Agreement validates this extractor stage only; it does not fix low GM/WM support or the fragmented output seen in the human quality audit.
- Distances use the user's symmetric nearest-surface limits of 0.01 mm mean and 1 mm Hausdorff; topology agreement is reported independently and is not inferred from distance agreement.
- The machine-readable full output, hashes, crop bounds, and runtime are in `build/invivo-boundary-jist-surface-comparison-2026-10-11/comparison.json` (ignored build output). The comparison used [`compare_in_vivo_jist_surfaces.py`](../../scripts/compare_in_vivo_jist_surfaces.py) with its core, surface, build, and report paths set to the corrected human tail locations.
