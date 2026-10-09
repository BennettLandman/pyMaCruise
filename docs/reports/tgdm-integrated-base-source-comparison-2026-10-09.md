# Integrated base-workflow TGDM core comparison against Java

## Result

The selected Java `MaGenericTGDM`/JIST calls were run on six held-fixed inputs
from the corrected seven-stage synthetic base workflow's `cruise-core`
stage. The Python comparator reran the selected TGDM core from those same
inputs and compared all eight saved stages against Java. All three delivered
surfaces pass the agreed criteria; tissue-sign output is exact.

| Stage | Zero-level Dice | Symmetric mean surface distance | Hausdorff | Result |
|---|---:|---:|---:|---|
| Central | 1.0 | `0.000001335 mm` | `0.003286 mm` | Pass |
| Inner | `0.99999894` | `0.000003170 mm` | `0.026402 mm` | Pass |
| Outer | `0.99998454` | `0.000161609 mm` | `0.088746 mm` | Pass |

The accepted limits are Dice ≥0.97, symmetric mean surface distance ≤0.01 mm,
and Hausdorff ≤1 mm. The Java/Python tissue-sign map is exact across all
5,317,312 voxels. Intermediate initial-inner, smoothed-initial, and
central-shrink sign maps have Dice ≥0.9999989 (the remaining differences are
recorded in the JSON report). The selected Python run reports all checks
passed.

This is a selected-source comparison with upstream inputs held fixed. It is
not a whole MaCRUISE pipeline comparison or real-MRI validation. The full
workflow outputs are in the original 172×216×164 grid; the TGDM reference
comparison is on the corresponding 166×208×154 RAS crop.

## Reproduction

```sh
BASE=build/base-workflow-2026-10-09-retained-center
./.venv/bin/python scripts/prepare_tgdm_core_java_inputs.py \
  "$BASE/stages/cruise-core/attempt-01" --target "$BASE/java-tgdm/input"
mkdir -p "$BASE/java-tgdm/reference" build/java-reference/tgdm-core-classes
javac -d build/java-reference/tgdm-core-classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/TGDMCoreReferenceHarness.java
/usr/bin/java -Xms512m -Xmx8g \
  -cp build/java-reference/tgdm-core-classes TGDMCoreReferenceHarness \
  "$BASE/java-tgdm/input" "$BASE/java-tgdm/reference"
./.venv/bin/python scripts/convert_tgdm_core_java_reference.py \
  --input-root "$BASE/java-tgdm/input" \
  --reference-root "$BASE/java-tgdm/reference"
./.venv/bin/python scripts/compare_tgdm_core_java_reference.py \
  --input-root "$BASE/java-tgdm/input" \
  --reference-root "$BASE/java-tgdm/reference" \
  --candidate-root "$BASE/java-tgdm/python-stage" \
  --report "$BASE/java-tgdm/comparison.json"
```

The exact metrics are also stored in
[`tgdm-integrated-base-source-comparison-2026-10-09.json`](tgdm-integrated-base-source-comparison-2026-10-09.json).
The Java harness completed successfully with an 8 GB maximum heap. The Python
comparison completed; its process RSS was not sampled during this invocation.
