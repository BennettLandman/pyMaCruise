# Matched Java/Python TOADS on human BASE inputs — 2026-10-13

## Result

Both naturally converged human TOADS passes match the selected Java wrapper
classification exactly when Java and Python receive the same working-space
intensity array and brain mask:

| Pass | Java / Python terminal iteration | Python terminal energy ratio | Changed labels / voxels | Per-label Dice |
|---|---:|---:|---:|---:|
| First TOADS | 21 / 21 | 4.86711651e-5 | 0 / 12,058,624 | 1.0 for all 11 labels |
| Second TOADS | 20 / 20 | 1.48473373e-5 | 0 / 12,058,624 | 1.0 for all 11 labels |

The outputs are byte-for-byte equal by voxel after reshaping Java's
x-fastest buffer in Fortran order. Each comparison covers labels
`1, 5, 6, 14, 15, 16, 17, 18, 22, 24, 25`; both input pairs have shape
`184×256×256`, 1 mm spacing, and RAS orientation. This is controlled
algorithm-level equivalence evidence for the two TOADS passes, not whole
MaCRUISE equivalence or independent anatomical acceptance.

The comparison deliberately separates TOADS from the historical human output
path. Java has no brain-mask argument in the selected wrapper harness, so the
same user mask is applied to the source N4 T1 before writing the Java
float32/x-fastest input; Python receives the original T1 and that mask and
applies it in its normal preparation path. Both then use the same bundled
VandyCRUISE descriptor and source-selected parameters. The second pair is the
Python Dura-stripped T1 plus its Dura brain mask, passed unchanged to Java and
Python.

## Held-fixed inputs and parameters

First pass:

- N4 input: `build/invivo-base-orientation-only-source-profile-2026-10-12/stages/n4/attempt-01/t1-n4.nii.gz` (SHA-256 `2236a4be0a25bbd4231994aad95fd614e7fee5043de289c5a3f2cc7c2f69cdae`)
- Brain mask: `build/invivo-base-orientation-only-source-profile-2026-10-12/resources/orientation-only-inputs/BrainMask.nii.gz` (SHA-256 `f45b1d4759098401eaddf65ffbd16f5a4dca051546ac1467af6c3608a3527672`)
- Applied-mask Java float32 input SHA-256: `31a350e070306aef86098704421366c55f0f2e7a4a15673e26774f98ce107638`
- Java classification SHA-256: `5c3ec8530dbeb2b222d7f35683c92e09e24dee6f12d2f92164c4cf8e0af5b165`
- Python NIfTI classification SHA-256: `611bcf86ac174b4cf133ea15b657e2ddb1548353c628e8e59be849f023478b27`

Second pass:

- Dura-stripped T1 SHA-256: `57f289a0e2916f79cc2b503fe9c0c0d29b0425e401c74e8435be0cf7b5581334`
- Dura brain mask SHA-256: `f4b50c67394f5a6956b30900ebe6442a974f2d6d8206e20b938a798d99e898e9` (1,277,381 positive voxels)
- Applied-mask Java float32 input SHA-256: `0ff0fc91c0ba4bc3d9abac5b811a0fb4fc04b0dd2188b559c491171668201bed`
- Java classification SHA-256: `18345f4a8a4cf29288cd9be2deb3f54b27e3170866d47d3977da5150e9348b90`
- Python NIfTI classification SHA-256: `585e960ebc016c0fbb2099c3a08a97ef306113401b96f9792cbca72037d917e6`

Both passes used `T1_MPRAGE`, 4 registration levels, 50 registration
iterations, 5 post-refresh registration iterations, registration maximum 5,
segmentation maximum 50, convergence threshold `1e-4`, maximum distance
`1e-4`, padding 10, and relation limit 3. The selected Java wrapper harness is
`tests/java_reference/ToadsWorkflowSourceHarness.java` (SHA-256
`397cb71fc3a8a18186fc132321bae54b7aa3691a9fd5331725f1b8287885c2c2`). The
atlas descriptor used by both implementations has SHA-256
`cfe3e2556cca554122f5f16427878b6924cb1de6947e9392b2918bed24a4a939`.

Python TOADS peak sampled RSS was 6.834 GB (first pass) and 7.043 GB (second
pass), each under a 16 GB watchdog. The selected Java harness stayed below
4.6 GB. The Python Dura-removal stage completed in between the two TOADS runs;
its stripped T1, brain mask, corrected GM/WM, and topology products reproduce
the prior full human Python BASE run byte-for-byte. The first- and second-pass
Python TOADS label maps also reproduce that run's maps byte-for-byte. This
confirms restart-stage repeatability for this code/data path; Dura itself was
not re-compared to Java in this experiment.

## Interpretation and remaining human BASE gaps

This resolves the prior uncertainty about the Python TOADS evolution on the
user-mask/native-RAS working inputs: the two TOADS passes match the selected
Java implementation at natural convergence, including every final label.
It does **not** resolve why the saved historical MaCRUISE output differs.
That run used a MaSeg-derived MultiAtlas mask and SPECTRE affine preprocessing
before TOADS, whereas this selected Python base workflow uses the supplied
T1-matched brain mask in the native T1 grid. The source output's missing
transform/provenance still prevents a held-fixed comparison through those
upstream stages.

The full Python BASE run already completed on this subject, but its remaining
acceptance gaps are downstream: refined-label Dice, surface-position error,
GM/WM ribbon support, and refinement provenance. This TOADS result does not
change the source-faithful refinement default, establish that the user mask
and no-regrid profile are the historical settings, or validate the meshes.
The next BASE work should resume at the downstream checkpoints and report
label support, surface topology/nesting, thickness, and any available source
comparison without rerunning the two now-matched TOADS passes.

## Reproduction

All generated Java raw arrays, NIfTIs, and manifests remain in ignored `build/`
under `human-toads-matched-2026-10-13`; no human images or outputs are
committed. The first and second comparison JSON files are stored in that same
directory. The Java wrapper harness and comparator are run with the selected
descriptor and source settings above. The run's Python command is recorded in
the previous full-run manifest at
`build/invivo-base-orientation-only-source-profile-2026-10-12/run.json`.

The successful result is source equivalence for selected TOADS on two matched
human inputs. It is neither synthetic evidence nor full-workflow anatomical
validation.
