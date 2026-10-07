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

## Current selected-core rerun

The earlier saved base-run comparison below recorded an outer Hausdorff failure
from the then-current Python implementation. That saved output predates the
source-semantics fix and is retained as historical evidence. The selected-core
comparison has now been rerun on the committed synthetic inputs by
`scripts/compare_tgdm_core_java_reference.py`; its current measurements are in
[the current stage JSON](tgdm-core-brain-resolution-java-comparison-2026-10-07.json).
The rerun passes: tissue-sign output is voxel-exact; inner/central/outer Dice
are `1.0/1.0/0.9999887`, symmetric mean surface distances are
`0.00000360/0.00000137/0.00014844 mm`, and Hausdorff distances are
`0.00243/0.00329/0.19498 mm`. The prior isolated outer discrepancy was due to
counting components not seeded from the center's 6/18-neighbor set. Java's
float32 subtraction order for the tissue-sign denominator also matters at two
quantization boundaries; Python now matches both behaviors. The isolated
point-level trace and focused TGDM tests are documented in
[the outer comparison report](tgdm-outer-shared-central-java-comparison-2026-10-07.md).

This direct selected-core result does not retroactively validate the old saved
base-run outputs or the full six-stage wrapper. The restartable base workflow
must be rerun with the correction before claiming updated end-to-end workflow
evidence.

The comparator was rerun on committed code `6e2e748`; its full JSON metrics are
in [the machine-readable report](tgdm-core-saved-base-run-java-comparison-2026-10-07.json).

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
