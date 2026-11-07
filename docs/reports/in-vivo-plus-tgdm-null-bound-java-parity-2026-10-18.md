# Human PLUS null-bound TGDM comparison — 2026-10-18

## Scope

This is a held-fixed comparison of the selected Java `MaGenericTGDM` sequence
and the Python null-bound TGDM route on case 100260's 140×178×134 RAS crop.
Both receive the same initial membership, tissue, ACE, skeleton, and GVF
arrays. Java performs the selected initial-membership conversion and signed
distance initialization. The optional `Brain Mask Level Set` is omitted from
the Java calls and from the Python TGDM call, matching the selected
`run_TGDMinvcrop_jist.m` wrapper. The supplied brain mask remains recorded as
an input/crop resource.

The Python input preparation records all six resulting stage arrays as exact
matches to the saved workflow outputs (`saved_stage_difference_voxels` is
empty). The held-fixed inputs are identified by SHA-256 in the
[machine-readable report](in-vivo-plus-tgdm-null-bound-java-parity-2026-10-18.json)
and the build-side `input-provenance.json`. The Java run used the previously
compiled reference class set; the main harness class hash is recorded in the
JSON. Java completed under the 8 GB watchdog at 1.259 GB peak sampled RSS.

## Results

The initial signed level set differs by at most 4.77e-7, with no sign changes.
All five zero-level surfaces passed the established held-fixed TGDM gates:

| Surface | Zero-level Dice | Symmetric mean distance (mm) | Hausdorff (mm) |
|---|---:|---:|---:|
| Smoothed initial | 1.000000 | 3.65e-7 | 0.006424 |
| Initial inner | 1.000000 | 2.18e-7 | 0.005065 |
| Central | 1.000000 | 1.51e-7 | 0.004322 |
| Inner | 1.000000 | 7.68e-8 | 0.000766 |
| Outer | 0.999994 | 1.10e-4 | 0.225497 |

Each surface is within Dice ≥0.97, symmetric mean distance ≤0.01 mm, and
Hausdorff ≤1 mm. Scalar fields are not bit-identical; their voxelwise errors
are reported in the JSON.

## Limits

This comparison holds the Python-generated upstream membership and fields
fixed. It does not validate AtlasEM, Dura, LesionTOADS, Aladin registration,
the selected atlas initialization, or an independently run complete PLUS
workflow. The result is TGDM-stage agreement for the selected null-bound
handoff only; whole-workflow PLUS equivalence remains unverified.
