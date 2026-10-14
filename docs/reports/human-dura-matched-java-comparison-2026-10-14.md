# Human Dura stage comparison against selected Java source (2026-10-14)

## Result

The Python Dura-removal stage matches the selected original Java Dura wrapper on the same held-fixed first-pass inputs from human BASE case 100260. This is a stage-level source comparison; it does not establish two-pass TOADS or whole-MaCRUISE equivalence.

The Java harness executes the selected `MedicAlgorithmRemoveDura09` sequence with the original topology-correction and TGDM implementations. Both implementations used the same 184 × 256 × 256 RAS inputs, five outer TGDM passes, and MaCRUISE wrapper membership level 0.3. The source harness supports alternate membership values so this run uses the wrapper-selected 0.3 rather than the synthetic harness's historical default 0.4.

## Input identity

The input NIfTIs are from `build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/` in the local validation checkout. They are held fixed between implementations and are not copied into this repository.

| Input | SHA-256 |
| --- | --- |
| T1 | `2236a4be0a25bbd4231994aad95fd614e7fee5043de289c5a3f2cc7c2f69cdae` |
| White matter | `f7d09ef1706669ee6b7d3d930f2efc6175342d8fb1d069f1bcf9a0703100f779` |
| Gray matter | `f8bdf9ce287ffdc670a46428f93a75152b57295a837a98b2ebe7889a29cc5b2f` |
| White-matter mask | `e26b8565bc8a2782cc874ae5b4ceb600b773ac48a6fe2739db8926e5eecd1034` |

## Comparison

The two topology-correction maps, their union, stripped T1, corrected WM/GM, and final brain mask are bit-identical. The initial inner and outer distance fields differ by at most `9.54e-7`; neither has a sign disagreement. All five TGDM pass fields retain identical nonpositive regions (Dice 1.0, zero sign disagreements); their largest absolute value difference is `4.34e-5`. The final Dura level set differs by at most `2.91e-5`, with zero sign disagreements, and the final brain-mask voxel arrays are identical.

For a geometry check, both final level sets were passed through the same Python source-derived surface extractor. Each mesh has 188,850 vertices and 377,696 faces. Symmetric vertex-to-triangle mean distance is `1.13e-8 mm`, and sampled Hausdorff distance is `3.44e-5 mm`. This isolates level-set-induced geometry under one extractor; it is not a comparison of Java `IsoSurfaceOnGrid` against the Python mesh extractor.

The Python comparison sampled peak RSS at 2.667 GB under a 12 GB watchdog. Focused Dura tests passed: **29 passed**, with two expected RAS reorientation warnings. The full suite was not rerun because this update changes only a validation script and documentation; the latest full-suite result for the unchanged production implementation is documented in [the 2026-10-14 refinement report](full-pytest-2026-10-14-refinement-human-refresh.md).

## Handoff into the second TOADS pass

The Java Dura run's `stripped_t1.raw.gz` and `brain_mask.raw.gz` were also
compared with the NIfTI handoff saved by the earlier Python BASE run. Both
arrays are voxel-identical (0 changed voxels each); the saved Python NIfTI
hashes are `57f289a0e2916f79cc2b503fe9c0c0d29b0425e401c74e8435be0cf7b5581334`
and `f4b50c67394f5a6956b30900ebe6442a974f2d6d8206e20b938a798d99e898e9`.
Those exact NIfTIs are the held-fixed inputs in the human second-pass
Java/Python comparison, where both natural-convergence classifications match
exactly. Together these results validate the Dura-to-second-TOADS handoff and
second-pass label result on this case. This composition still begins from
held-fixed first-pass tissue maps; it does not validate first-pass TOADS,
their upstream preprocessing, or the whole BASE workflow as one independent
Java execution.

## Reproduction

From the repository root, construct the Java input pack from the four files above, then run the selected Java harness:

```sh
DURA_WRAPPER_NX=184 DURA_WRAPPER_NY=256 DURA_WRAPPER_NZ=256 \
DURA_WRAPPER_ITERATIONS=5 DURA_WRAPPER_MEMBERSHIP=0.3 DURA_WRAPPER_JAVA_HEAP=8g \
  sh scripts/run_java_dura_base_reference.sh \
  build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-input.bin \
  build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-output-membership-03
```

Then compare the Java output with the Python stage on the identical inputs:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/compare_dura_java_outputs.py \
  --t1 build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/t1.nii.gz \
  --white-matter build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/white_matter.nii.gz \
  --gray-matter build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/gray_matter.nii.gz \
  --white-matter-mask build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/white_matter_mask.nii.gz \
  --java-output-dir build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-output-membership-03 \
  --membership-level 0.3 --iterations 5 \
  --report build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-comparison-membership-03.json
```

The generated input pack, Java outputs, JSON comparison, and human image files reside under ignored `build/` data and are not committed. The Java harness and reusable Python comparison command are committed with this report.

## Scope and open work

This validates one Dura invocation on the first-pass handoff. It does not validate SPECTRE registration, TOADS classification, the second TOADS pass, refinement, surfaces from the selected Java meshing implementation, or whole-workflow anatomy and label acceptance. The parent P18/PyMaCRUISE workflow remains open.
