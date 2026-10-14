# Human Dura stage comparison against selected Java source (2026-10-14)

## Result

The Python Dura-removal stage matches the selected Java call sequence on the same held-fixed first-pass inputs from human BASE case 100260. This comparison uses the MaCRUISE wrapper's actual topology regularization settings: WM `0.05`, GM `0.02`. It is stage-level source evidence, not whole-BASE equivalence.

The selected MATLAB wrapper is `pipeline/run_RemoveDura_jist.m` (SHA-256 `57473ab813e81b89ffed433de89b2fad8e111f6b68fa41e30c920639cd864e80`); it passes `-inWM3 0.05` and `-inGM2 0.02`. The corresponding `toads-cruise` `MedicAlgorithmRemoveDura09.java` exposes and applies those two values (SHA-256 `605c3ef4022d0eaabe8a609a6c08400dd8923c6075f9b4c4fbc34b500c114147`). The older VandyCRUISE copy (SHA-256 `1ed180daf4a4334c995ad310828ec91cfd5593413eec17ab15404a421d7805a4`) hardcodes 0.05 for both topology calls; it is not the MaCRUISE V3 parameter contract. Earlier Dura reports and references that used 0.05 for both are superseded for the V3 wrapper profile.

The Java harness executes the selected Dura call sequence with the original topology-correction, JIST distance-field, and TGDM implementations. It is not a launch of the full MIPAV/JIST wrapper application. Both implementations used the same 184 × 256 × 256 RAS inputs, five TGDM outer passes, and wrapper-selected membership level 0.3.

## Input identity

The held-fixed NIfTIs are in ignored `build/` output from the local validation checkout and are not committed.

| Input | SHA-256 |
| --- | --- |
| T1 | `2236a4be0a25bbd4231994aad95fd614e7fee5043de289c5a3f2cc7c2f69cdae` |
| White matter | `f7d09ef1706669ee6b7d3d930f2efc6175342d8fb1d069f1bcf9a0703100f779` |
| Gray matter | `f8bdf9ce287ffdc670a46428f93a75152b57295a837a98b2ebe7889a29cc5b2f` |
| White-matter mask | `e26b8565bc8a2782cc874ae5b4ceb600b773ac48a6fe2739db8926e5eecd1034` |

## Comparison

Both topology-correction outputs and their union are bit-identical. The inner and outer distance fields differ by at most `9.54e-7`; neither has a sign disagreement. Every one of the five TGDM pass fields has the same nonpositive region (Dice 1.0, zero sign disagreements); the largest pass-field absolute difference is `6.66e-5`. The final level set differs by at most `2.48e-5`, with zero sign disagreements. The final brain-mask arrays are bit-identical.

The final Java and Python level sets were passed through the same Python source-derived extractor. Both meshes have 175,384 vertices and 350,764 faces; symmetric vertex-to-triangle mean distance is `7.72e-9 mm`, and sampled Hausdorff distance is `3.05e-5 mm`. This checks geometry induced by the held-fixed fields under one extractor; it does not compare Java `IsoSurfaceOnGrid` to the Python extractor.

The Java output's stripped T1 and brain mask each differ at **89,885 voxels** from the previously saved Python handoff that was used in the older human second-TOADS comparison. Therefore that prior second-pass comparison is not a downstream validation of this corrected Dura profile. The second TOADS pass and the full no-regrid BASE run need refreshing with this handoff before claiming a composed human result.

The Python comparison sampled peak RSS at **2.721 GB** under a 12 GB watchdog. The focused Dura suite passed **32 tests**, with two expected RAS reorientation warnings.

## Reproduction

Build the Java input pack from the four held-fixed files, then run the selected Java harness with the wrapper values:

```sh
DURA_WRAPPER_NX=184 DURA_WRAPPER_NY=256 DURA_WRAPPER_NZ=256 \
DURA_WRAPPER_ITERATIONS=5 DURA_WRAPPER_MEMBERSHIP=0.3 \
DURA_WRAPPER_GM_REGULARIZE=0.02 DURA_WRAPPER_JAVA_HEAP=8g \
  sh scripts/run_java_dura_base_reference.sh \
  build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-input.bin \
  build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-output-membership-03-gm002
```

Compare selected Java and Python on the identical inputs and save the Python NIfTI products for the next stage:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/compare_dura_java_outputs.py \
  --t1 build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/t1.nii.gz \
  --white-matter build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/white_matter.nii.gz \
  --gray-matter build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/gray_matter.nii.gz \
  --white-matter-mask build/invivo-base-orientation-only-source-profile-2026-10-12/dura-java-first-pass-input/white_matter_mask.nii.gz \
  --java-output-dir build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-output-membership-03-gm002 \
  --python-output-dir build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-python-output-membership-03-gm002 \
  --membership-level 0.3 --iterations 5 \
  --report build/invivo-base-orientation-only-source-profile-2026-10-14/dura-human-java-comparison-membership-03-gm002.json
```

The input pack, subject images, Java/Python products, and JSON comparison remain under ignored `build/` data.

## Scope and open work

This validates Dura on held-fixed first-pass tissue maps. It does not validate SPECTRE registration, upstream TOADS, the second TOADS pass with the corrected 0.02-GM-regularization handoff, refinement, or whole-workflow anatomy and label acceptance. The parent P18/PyMaCRUISE workflow remains open.
