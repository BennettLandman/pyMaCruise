# Numerical Validation Tolerances

**Status:** P08 started. The project owner supplied the target tolerances on 2026-09-27. Metric definitions below make those targets operational; the surface-distance interpretation remains provisional pending owner confirmation and synthetic validation.

## Acceptance thresholds

| Output | Metric | Acceptance threshold | Units / interpretation |
|---|---|---|---|
| Segmentation masks | Dice similarity coefficient | At least `0.97` (allowing a 0.03 absolute Dice difference) for each required segmentation label/region. | Unitless; compare corresponding masks after confirming the same physical grid. Both-empty masks score 1; one-empty masks score 0. |
| Surfaces | Symmetric Hausdorff distance | At most `1.0` | Millimetres in RAS+; maximum of the two directed surface distances. |
| Surface positions | Mean absolute surface distance | At most `0.01` | Millimetres in RAS+; provisional definition is the mean of both directed nearest-point-to-surface distances, using physical surface geometry rather than raw voxel indices. |

## Comparison rules

- Express all spatial comparisons in RAS+ millimetres. Normalize non-RAS inputs before comparison and preserve the warning/conversion record.
- Compare each required segmentation label or region separately and report the per-label Dice values; also report the aggregate summary, without allowing a strong label to hide a weak one.
- Compute surface distances symmetrically so neither reference nor candidate is treated as the sole source of distance samples. Use a consistent surface-distance implementation and document its sampling or point-to-triangle method.
- Report the compared image grids, surface frames, transforms, parameters, and source/software versions with each validation result.
- Keep exact equality as the expected criterion for deterministic label operations whose source contract is exact. The Dice threshold applies to end-stage segmentation equivalence and other explicitly designated approximate comparisons.

## Pending confirmation and evidence

- Confirm that the 3% segmentation tolerance means Dice `>= 0.97` for each required label/region.
- Confirm that “mean absolute distance for surface positions” means the symmetric mean nearest-surface distance defined above, rather than pointwise vertex-coordinate distance.
- Validate the metric implementation on synthetic masks and known-position surfaces before using these thresholds for source equivalence or real-data acceptance.

These values are project acceptance targets, not claims that any current implementation has passed them. Stage-specific tolerances may be tighter when exact or deterministic behavior is expected; any justified exception must be recorded with evidence.
