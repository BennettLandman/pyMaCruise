# Numerical Validation Tolerances

**Status:** P08 complete. The project owner confirmed both per-region Dice ≥0.97 and symmetric mean nearest-surface distance for the 0.01 mm surface-position threshold on 2026-09-27. Synthetic comparator fixtures pass.

## Port acceptance objective — owner clarification, 2026-10-14

The objective is to reproduce the selected original Java/MATLAB MaCRUISE workflows within the agreed numerical tolerances, including their limitations. Improving anatomical quality beyond the original is outside the current port scope. Successful execution alone or selected-stage parity alone does not establish whole-workflow equivalence.

The thresholds below measure **Python-versus-reference agreement** on matched inputs, preprocessing, grids, atlas, parameters, source versions, and stopping rules. In particular, per-label Dice ≥0.97 is an equivalence requirement; it is not a requirement that 97% of supplied tissue labels lie within the generated surfaces.

Ribbon tissue support, connected components, Euler characteristic, nesting, and absolute thickness bounds remain diagnostic measurements. Existing audit `passed` fields refer to their stated diagnostic criteria, not overall port acceptance. A source-matched result may pass port acceptance even if both implementations fail those anatomical criteria. An anatomical failure by Python alone neither proves a port defect nor establishes parity: compare the original on the same handoff or complete workflow before assigning a cause. Known source behavior, including defects and iteration-cap stops, is retained and documented unless a separate change is explicitly requested.

Analytic synthetic fixtures still check the contracts they were designed to test. A fixture's prescribed thickness or topology must not become a universal human-workflow requirement. For human thickness, compare Python with the original measurement method and output using stage-specific tolerances; do not require a generic 5 mm thickness.

Prioritize reproducible stage and complete-workflow reference comparisons and investigate differences. If the reference or its settings are unavailable, report equivalence as unverified rather than substitute ideal anatomy for reference evidence. Historical audit measurements remain valid observations; their earlier wording about anatomical acceptance is superseded by this scope clarification.

## Acceptance thresholds

| Output | Metric | Acceptance threshold | Units / interpretation |
|---|---|---|---|
| Segmentation masks | Dice similarity coefficient | At least `0.97` for each required segmentation label/region. | Unitless; compare corresponding masks after confirming the same physical grid. Both-empty masks score 1; one-empty masks score 0. |
| Surfaces | Symmetric Hausdorff distance | At most `1.0` | Millimetres in RAS+; maximum of the two directed surface distances. |
| Surface positions | Symmetric mean nearest-surface distance | At most `0.01` | Millimetres in RAS+; mean of both directed vertex-sample-to-triangle-surface distances. |

## Comparison rules

- Express all spatial comparisons in RAS+ millimetres. Normalize non-RAS inputs before comparison and preserve the warning/conversion record.
- Compare each required segmentation label or region separately and report the per-label Dice values; also report the aggregate summary, without allowing a strong label to hide a weak one.
- Compute surface distances symmetrically so neither reference nor candidate is treated as the sole source of distance samples. The implementation computes exact closest points on triangle interiors/edges and uses each mesh's vertices as deterministic surface samples.
- Report the compared image grids, surface frames, transforms, parameters, and source/software versions with each validation result.
- Keep exact equality as the expected criterion for deterministic label operations whose source contract is exact. The Dice threshold applies to end-stage segmentation equivalence and other explicitly designated approximate comparisons.

## Output-class acceptance methods

| Output class | Primary comparison | Supporting checks |
|---|---|---|
| Discrete labels and binary masks | Exact voxel equality for deterministic stages; per-label Dice ≥0.97 for final approximate equivalence. | Label set, voxel counts, grid/affine equality, and RAS axis codes. |
| Continuous scalar volumes | MAE, RMSE, maximum absolute error, plus exact checks where the stage is deterministic. | Grid/affine, dtype, finite-value policy, and units. P08 establishes no universal scalar threshold; set stage-specific limits before equivalence claims. |
| Surfaces and point sets | Symmetric triangle-surface Hausdorff ≤1 mm and mean nearest-surface distance ≤0.01 mm in RAS+ mm. | Watertightness, connected components, Euler characteristic, bounds, normals, and units. |
| Rigid/affine transforms | Translation error in mm and rotation error in degrees. | Matrix validity, source/target frame, direction, and resampling reference grid. Stage-specific numeric limits are set before equivalence claims. |
| Thickness and derived scalar measurements | Voxel/vertex MAE, RMSE, maximum error, and distribution summaries. | Valid/invalid counts, units, mesh/voxel frame, and stage-specific limits set before equivalence claims. |

## Confirmation and evidence

- The owner confirmed that the 3% segmentation tolerance means Dice `>= 0.97` for every required label/region.
- The owner confirmed the symmetric mean nearest-surface interpretation. Synthetic tests verify zero distance for identical surfaces with different triangulations and the expected metrics for a known 0.02 mm offset.
- `tests/test_validation.py` exercises empty/nonempty and per-label Dice, deliberately perturbed voxels, symmetric point and triangle-surface distances, transform error metrics, topology checks, and deterministic JSON report output.

These values are project acceptance targets, not claims that any current implementation has passed them. Stage-specific tolerances may be tighter when exact or deterministic behavior is expected; any justified exception must be recorded with evidence.
