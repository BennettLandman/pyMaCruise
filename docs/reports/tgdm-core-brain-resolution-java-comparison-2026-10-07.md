# CRUISE core TGDM Java comparison — 2026-10-07

## Scope and command

This is a selected-source comparison on the committed synthetic atlas-derived
brain-resolution fixture. It contains no human MRI. The Java harness executes
the selected `MaGenericTGDM` and `DistanceField` calls for the same crop inputs,
nonzero GVF, and one-iteration wrapper profile recorded in the fixture
manifest. The saved Python results come from the restartable synthetic base
workflow's `toads-cruise-core` stage. Java crop outputs are transferred back to
the native RAS grid with the selected JIST-style legacy uncrop rule before
comparing the three delivered cortical surfaces.

Fixture and Java reference regeneration:

```bash
./.venv/bin/python scripts/prepare_tgdm_core_java_inputs.py \
  build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-02
sh scripts/generate_java_tgdm_core_reference.sh
```

Comparison command:

```bash
./.venv/bin/python scripts/compare_tgdm_core_saved_java_reference.py \
  --report docs/reports/tgdm-core-saved-base-run-java-comparison-2026-10-07.json
```

The report command exits nonzero because the outer-surface Hausdorff gate is
missed. This is a concrete integration blocker, not evidence that the complete
workflow is source-equivalent.

## Results

| Output | Dice | Symmetric mean surface distance (mm) | Hausdorff (mm) | Gate |
|---|---:|---:|---:|---|
| Inner | 0.9999634 | 0.0003150 | 0.763979 | Pass |
| Central | 0.9999980 | 0.0000103 | 0.300246 | Pass |
| Outer | 0.9999654 | 0.0005205 | 1.562443 | **Fail: >1 mm** |

The outer masks differ at 43 of 6,092,928 native-grid voxels. Eight level-set
voxels differ by more than 1.0; the maximum absolute level-set difference is
1.7052. The Hausdorff outlier is localized near native voxel `(95, 87, 66)`.
The outer output therefore does not meet the agreed 1 mm maximum, despite its
high Dice and small mean surface distance. Do not treat this as a passed
CRUISE core equivalence check.

## Interpretation and next work

The Java source call sequence, input affine, crop bounds, RAS orientation,
vector-field dimensions, and selected one-iteration settings are captured in
the committed harness and manifests. Inner and central pass the agreed gates.
The single localized outer-surface deviation is now the next TGDM integration
investigation: compare the Java/Python central level sets, outer force inputs,
and final narrow-band updates around the reported voxel. Keep P34/P37 Started
until this discrepancy is explained and the required workflow-level checks
pass. No real-data validation was performed.
