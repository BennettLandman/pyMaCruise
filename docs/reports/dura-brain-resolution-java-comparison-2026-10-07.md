# Brain-resolution base-run Dura stage comparison — 2026-10-07

## Scope

This compares the selected `MedicAlgorithmRemoveDura09` wrapper sequence on
the exact synthetic maps fed to the corrected brain-resolution base workflow:
N4-corrected T1 plus first-pass TOADS WM-fill, GM, and WM mask. The maps are
held fixed for Java and Python. They are RAS NIfTIs on the `(172,216,164)`
grid at 0.9375 mm; no human scan is used. The inputs and selected Java outputs
are committed under
[`tests/data/synthetic/base_workflow_brain_resolution/dura_java_inputs`](../../tests/data/synthetic/base_workflow_brain_resolution/dura_java_inputs)
with hashes and provenance.

The selected Java wrapper sequence runs both 6/18 topology corrections, the
WM+GM union, initial inner/outer distance fields, and one Dura TGDM outer pass
(30 internal updates, matching the current bounded base-run setting), then
thresholds the mask at 2 and masks the returned T1/WM/GM. The Java harness
uses the copied selected `AlgorithmTopologyCorrection`, `DistanceFieldyh`
(the supplied source states it is identical to the omitted `DistanceField`),
and selected VandyCRUISE TGDM sources. It is a source-based headless harness,
not a full MIPAV/JIST application run.

## Results

The two topology-corrected maps and their WM+GM union match exactly. The
initial inner and outer distance fields have zero sign disagreements; their
maximum absolute differences are `0.01394` and `0.01122` voxels, while their
mean absolute errors are `1.22e-8` and `1.84e-8`.

The one-pass TGDM level set differs at 323,396 voxels, with maximum absolute
difference `0.55564` and mean `4.63e-5`. There are 18 zero-level-set sign
disagreements out of 6,092,928 voxels; Dice is **0.999990**. At the wrapper's
actual mask threshold of 2, the physical surface has symmetric mean distance
**0.000308 mm** and Hausdorff distance **0.159904 mm**.

The final brain mask differs at **23 of 6,092,928 voxels**. Foreground Dice is
**0.9999896** (1,101,377 Java voxels; 1,101,356 Python voxels). Its physical
surface has symmetric mean distance **0.000214 mm** and Hausdorff distance
**0.9375 mm**, within the agreed Dice ≥0.97, mean distance ≤0.01 mm, and
Hausdorff ≤1 mm gates. The same 23 voxels cause differences in masked T1 and
WM/GM outputs. The GM membership has a large per-voxel maximum difference
where one output is masked and the other is not; the difference is localized
to those 23 mask-boundary voxels, not distributed through the tissue field.

The Python comparator also confirms its regenerated outputs are byte-for-byte
identical to the six-stage synthetic base run's saved Dura outputs. Thus this
run both reproduces the existing stage and compares that stage to selected
Java on the same inputs.

## Reproduction

```sh
./.venv/bin/python scripts/prepare_dura_base_java_input.py \
  tests/data/synthetic/base_workflow_brain_resolution/dura_java_inputs \
  build/java-reference/dura-base-input.bin

./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  scripts/run_java_dura_base_reference.sh \
  build/java-reference/dura-base-input.bin \
  build/java-reference/dura-base-output

./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_dura_base_java.py \
  build/java-reference/dura-base-output \
  --base-stage build/ace-indexed-heap-2026-10-07/base-run/stages/dura/attempt-02
```

The Java runner has an 8 GB maximum heap and ran the selected topology and
single-pass TGDM sequence to completion. The Python comparison peaked at
**2.107 GB** under a 16 GB watchdog. The committed `.nii.gz` references are
the actual Java arrays; temporary raw/gzip interchange files remain under the
ignored `build/` directory.

## Boundary and next integration work

This validates the Dura stage on the base workflow's actual synthetic inputs;
it does not independently validate the upstream TOADS outputs that produced
those inputs, the second TOADS pass, or the whole MaCRUISE pipeline. The 23
mask-boundary voxel discrepancy is retained explicitly. The Hausdorff result
is within 1 mm but relatively close to that limit, so retain it as a regression
metric. Next compare the complete base chain at the CRUISE core and final
refinement outputs, then use the approved real-data run only after the
synthetic chain is source-consistent.
