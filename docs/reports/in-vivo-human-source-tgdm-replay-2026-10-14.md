# Downloaded human BASE TGDM source replay

**Evidence category:** held-fixed Java source-stage replay on the downloaded human MaCRUISE outputs. This is not a fresh full Java/MATLAB workflow run.

## Result

The selected Java TGDM sequence reproduces the downloaded inner, central, and outer level-set surfaces exactly. For each output, Java and the downloaded source output have identical zero-level sign masks, Dice **1.0**, symmetric mean surface distance **0 mm**, sampled Hausdorff distance **0 mm**, and identical extracted vertex and face counts. The full float fields differ away from the surface, but no voxel sign changes occur:

| Output | Different field voxels | Mean absolute field error | Maximum absolute field error | Sign differences | Dice | Mean surface distance | Hausdorff |
|---|---:|---:|---:|---:|---:|---:|---:|
| Central | 19,747 | 0.0371 | 14.5947 | 0 | 1.0 | 0 mm | 0 mm |
| Inner | 16,823 | 0.0275 | 14.0210 | 0 | 1.0 | 0 mm | 0 mm |
| Outer | 24,699 | 0.0569 | 15.5974 | 0 | 1.0 | 0 mm | 0 mm |

Distances use the source Raw/T1 affine to express cropped-grid surfaces in RAS millimeters. Crop bounds were recovered from the downloaded `CorticalSeperation/Atlas_GM.nii.gz` using the selected 5-voxel padding rule: inclusive bounds `(21–160, 33–201, 71–204)`, shape `140×169×134`. The six exported TGDM inputs were all on the same full source grid and are hash-recorded in the [machine-readable report](in-vivo-human-source-tgdm-replay-2026-10-14.json).

## Reconstructing the missing GVF vector

The downloaded file named `TGDMinputs/target_image_GVFimg_uc_invcp_refspace.nii.gz` is scalar, while the selected TGDM wrapper consumes a three-component GVF. The replay therefore regenerated a vector from the exported ACE-enhanced-GM crop using selected Java `GradVecFlowGeneralized` with lambda `1.0`, 200 iterations, and normalization off. The exported scalar is most consistent with component 0 (correlation `0.9874`, MAE `0.00412`); components 1 and 2 have correlations near zero. The scalar export is insufficient to recover the original vector exactly, so this result does not claim full TGDM input identity. Even with this reconstructed vector, all three resulting source surfaces match the downloaded surfaces exactly.

The harness compiles the committed source copies of `MaGenericTGDM`, `DistanceFieldyh`, and `GradVecFlowGeneralized`; each copy's SHA-256 is checked against the supplied source checkout and recorded in JSON. The run uses the selected `run_TGDMinvcrop_jist.m` parameter values and call sequence, but does not execute the MIPAV/JIST CLI wrapper class itself. Whole-workflow equivalence, including the upstream TOADS, ACE, inputs, and exact historical GVF, remains open.

## Python versus Java on identical inputs

The replay now runs `tgdm_nested_arrays` on the same six cropped scalar inputs and regenerated three-component Java GVF used by the Java harness, with the selected WM iso level `0.5`. This removes NIfTI adaptation, vector regeneration, crop construction, and upstream workflow differences from the Python-versus-Java TGDM comparison.

| Output | Differing field voxels | Mean absolute field error | Maximum field error | Sign differences | Dice | Symmetric mean surface distance | Sampled Hausdorff |
|---|---:|---:|---:|---:|---:|---:|---:|
| Central | 87,418 | 3.51×10⁻⁹ | 1.91×10⁻⁶ | 0 | 1.000000 | 3.43×10⁻⁹ mm | 0.000016 mm |
| Inner | 92,274 | 4.72×10⁻⁹ | 4.50×10⁻⁶ | 0 | 1.000000 | 3.91×10⁻⁹ mm | 0.000015 mm |
| Outer | 2,179,517 | 3.08×10⁻⁴ | 0.4471 | 33 | 0.999984 | 0.000324 mm | 0.4541 mm |

All three surfaces meet the agreed Dice ≥0.97, mean surface distance ≤0.01 mm, and Hausdorff ≤1 mm limits. Inner and central level sets differ only at float precision; outer has a small set of sign changes and still passes the surface-position gates. The JSON report includes exact counts, mesh sizes, and metrics. This is held-fixed selected-source stage parity on human-derived inputs; the reconstructed historical GVF and missing independent original full workflow still limit whole-run conclusions.

The complete replay, including the Python evolution and mesh comparisons, peaked at **2.831 GB** under an **8 GB** RSS watchdog. This is strong human-data evidence for the selected TGDM surface stage, not clinical validation or full MaCRUISE acceptance.

## Reproduction

With the supplied source checkout and downloaded output tree available:

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/replay_human_source_tgdm.py \
  --work-dir build/human-source-tgdm-replay-2026-10-14 \
  --report docs/reports/in-vivo-human-source-tgdm-replay-2026-10-14.json
```
