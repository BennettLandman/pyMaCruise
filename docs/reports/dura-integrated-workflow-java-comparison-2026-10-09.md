# Integrated base-workflow Dura stage comparison against Java

## Result

The selected Java Dura wrapper was rerun on the exact four inputs emitted by
the corrected synthetic base workflow's first TOADS stage. This compares a
held-fixed Dura stage from the current integration, not upstream or full
workflow equivalence.

| Output | Java/Python comparison |
|---|---|
| White-matter topology | Exact, 0 / 6,092,928 changed voxels |
| Gray-matter topology | Exact, 0 / 6,092,928 changed voxels |
| Brain mask | 3 differing voxels; foreground Dice `0.9999986109` |
| Brain-mask level-set sign | 4 sign-different voxels; MAE `1.31677e-5`, maximum absolute error `1.38344` |
| Stripped T1 | 3 differing voxels; MAE `0.000184314` |
| Corrected white matter | 3 differing voxels; MAE `1.30992e-7` |
| Corrected gray matter | 3 differing voxels; MAE `0.000125408` |

The three differing mask voxels are Java-positive/Python-zero at `(42,85,110)`,
`(47,68,29)`, and `(62,62,114)` in NIfTI array order. A symmetric distance
transform comparison of 6-connected boundary-voxel-center samples gives mean
nearest-boundary distance `0.0000443541 mm` and maximum `0.9375 mm`. This is a
discrete voxel-boundary sample metric, not a continuous triangle-mesh
measurement. All results are on a synthetic RAS grid with isotropic 0.9375 mm
voxels. This stage meets the agreed Dice threshold; the four level-set sign
differences and three output-voxel differences remain explicit residuals.

| Provenance | Value |
|---|---|
| Python implementation | `b785bff` (`Match retained TOADS registration center`) |
| Shape / spacing | `(172, 216, 164)` / `0.9375 mm` isotropic |
| Java input SHA-256 | `a79f0dbe332b3e61296dee8b0c76016a87389f0587505afc45a6a13f3b5cd5cc` |
| Java brain-mask output SHA-256 | `a06a3a8a12b62a2cc5e27e0f09cdb0c8ccdf04932f4c429e4d442563f23b3db7` |
| Python brain-mask SHA-256 | `0421233b26d62ff540755949a19768a711c235d01ce618f80b3dba8273b63bfb` |
| Java wrapper | `DuraBaseInputsReferenceHarness` / `AlgorithmTopologyCorrection` |
| Validation | Java harness exited 0 under a 16 GB watchdog; one selected TGDM pass completed |

## Reproduction

The binary input has a big-endian three-int dimension header followed by four
big-endian float32 arrays in C order (T1, white-matter fill, gray matter,
white-matter mask). The four source images are the saved current-run products
under `build/base-workflow-2026-10-09-retained-center/stages/n4/attempt-01/`
and `stages/toads-first/attempt-01/dura-inputs/`.

```sh
./.venv/bin/python scripts/pack_dura_base_java_inputs.py \
  build/base-workflow-2026-10-09-retained-center/stages/n4/attempt-01/t1-n4.nii.gz \
  build/base-workflow-2026-10-09-retained-center/stages/toads-first/attempt-01/dura-inputs/t1-n4_wmfill.nii.gz \
  build/base-workflow-2026-10-09-retained-center/stages/toads-first/attempt-01/dura-inputs/t1-n4_gm.nii.gz \
  build/base-workflow-2026-10-09-retained-center/stages/toads-first/attempt-01/dura-inputs/t1-n4_wmmask.nii.gz \
  build/base-workflow-2026-10-09-retained-center/java-dura/input.bin

./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  env DURA_WRAPPER_JAVA_HEAP=8g \
  scripts/run_java_dura_base_reference.sh \
  build/base-workflow-2026-10-09-retained-center/java-dura/input.bin \
  build/base-workflow-2026-10-09-retained-center/java-dura/output

./.venv/bin/python scripts/compare_dura_current_stage.py \
  build/base-workflow-2026-10-09-retained-center/java-dura/output \
  build/base-workflow-2026-10-09-retained-center/stages/dura/attempt-01
```

The source inputs are preserved in the completed run manifest. Java/Python
arrays are compared with the big-endian gzip reader in
[`compare_dura_current_stage.py`](../../scripts/compare_dura_current_stage.py).
The comparator's committed-fixture mode in
[`compare_dura_base_java.py`](../../scripts/compare_dura_base_java.py) is
deliberately separate and should not be mistaken for this integrated-input
replay.
