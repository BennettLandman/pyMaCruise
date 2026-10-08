# Current seven-stage base run: selected Java TGDM comparison

## Scope

This report compares the CRUISE core outputs from the complete seven-stage
synthetic base workflow run recorded in
[`base-workflow-seven-stage-brain-resolution-2026-10-08.md`](base-workflow-seven-stage-brain-resolution-2026-10-08.md)
against a fresh run of the selected Java `MaGenericTGDM` harness using the
exact six saved inputs from that workflow. This validates a selected synthetic
source slice on the current workflow output. It does not establish whole
MATLAB/JIST pipeline equivalence, an independent MATLAB MaSeg-refinement
comparison, or real-MRI behavior.

The Java harness output arrays, after uncropping to the native grid, were
bit-for-bit identical to the previously committed selected-Java reference
arrays. The current run's core inputs also match the committed input manifest
by SHA-256. The Python CRUISE outputs themselves are separate current-run
files, identified below.

## Comparison results

| Output | Dice | Symmetric mean surface distance (mm) | Hausdorff (mm) |
|---|---:|---:|---:|
| Inner | 1.000000000 | 0.000003533 | 0.002414 |
| Central | 1.000000000 | 0.000001369 | 0.003286 |
| Outer | 0.999989543 | 0.000147829 | 0.200093 |

All three selected-source outputs pass the agreed gates: Dice at least 0.97,
symmetric mean nearest-surface distance at most 0.01 mm, and Hausdorff distance
at most 1 mm. The comparison used the native `(172,216,164)` RAS grid at
0.9375 mm spacing. The outer scalar field is not voxel-identical; the surface
metrics describe the resulting boundary agreement.

The source-faithful Python MaSeg-refinement rule applied to the fresh Java
inner/outer fields differs from the workflow's saved `SEG_refine.nii.gz` at
seven voxels. Labels 1 and 4 are exact; label 100 has Dice `0.9999751513`.
Those seven output differences are exactly the original label-100 voxels where
the Java and Python outer fields fall on opposite sides of the source
`0.2 mm / 0.9375 mm = 0.2133333333` voxel threshold. This attributes the label
residual to the upstream outer-field difference propagating through the same
source-faithful Python refinement operation. It is not an independent Java or
MATLAB refinement-output comparison.

## Inputs and provenance

Committed synthetic input hashes:

| Input | SHA-256 |
|---|---|
| T1 | `ebf945d06cc2519114763c949af0fc1ec88d938677748a41126852ab6782f158` |
| Brain mask | `2d617064fbea8e8a5d807165d2567d87a9b40f0d4487078810b25a755ad87a00` |
| MaSeg | `06cb8f49e12d7966c6ff025601695b41de9dae5526358ab23b2be7cb4a2631005` |

Current run outputs (under ignored `build/`):

| Output | SHA-256 |
|---|---|
| Python inner | `31233e72558c488a7e536373571c12293f9ca64d464e84f324c0f4cff056d042` |
| Python central | `4fdaf13b1d1ffd80df8e7c0e2a6c3b2290cc761da16148987166ecfa6aa0445c` |
| Python outer | `1f88ba7b4483fa6e59af759055db54577046cefcb5516045e5b7d469a3ace6c2` |
| Python refined MaSeg | `689a31da970b635699888452bff9043283c682daf580483c3662ae2a3801bbc6` |

Harness source hashes: `TGDMCoreReferenceHarness.java` is
`04284ad617d933f61844a222ea2e45dc0121b4ec94160c8b0e8d3342e05cd3c4`; selected
`MaGenericTGDM.java` is
`66124d4180864edc77f6fff280ab817eb0fffc061fdf427c9d8d725f850b211d`, and
`GenericTGDM.java` is
`70f5b21541e8a6c6139f1bd227c4ece12e8818db717f2a78ef50d69fe89f08d1`.
The input manifest and Java outputs were generated under ignored `build/`;
they need not be committed because the synthetic NIfTI inputs, selected-source
harness, and matching prior Java reference arrays are already committed.

## Reproduction

From the repository root, first run the seven-stage command in the linked
workflow report. Then:

```sh
./.venv/bin/python scripts/prepare_tgdm_core_java_inputs.py \
  build/ace-indexed-heap-2026-10-08/base-seven-stage-run/stages/cruise-core/attempt-01 \
  --target build/tgdm-current-seven-stage/java_inputs
mkdir -p build/tgdm-current-seven-stage/classes
javac -d build/tgdm-current-seven-stage/classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/TGDMCoreReferenceHarness.java
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 0.5 --report-interval-seconds 10 -- \
  java -Xms256m -Xmx3g -cp build/tgdm-current-seven-stage/classes \
  TGDMCoreReferenceHarness build/tgdm-current-seven-stage/java_inputs \
  build/tgdm-current-seven-stage/java_reference
./.venv/bin/python scripts/convert_tgdm_core_java_reference.py \
  --input-root build/tgdm-current-seven-stage/java_inputs \
  --reference-root build/tgdm-current-seven-stage/java_reference
./.venv/bin/python scripts/compare_tgdm_core_saved_java_reference.py \
  --core-dir build/ace-indexed-heap-2026-10-08/base-seven-stage-run/stages/cruise-core/attempt-01 \
  --reference-root build/tgdm-current-seven-stage/java_reference \
  --report build/tgdm-current-seven-stage/comparison.json
```

The Java harness exited 0 and peaked at **1.605 GB sampled RSS** under its
4 GB watchdog. The comparator exited 0 with `passed=true`. This source
comparison did not run pytest or the full suite; the suite evidence for the
implementation is recorded in the linked workflow report.
