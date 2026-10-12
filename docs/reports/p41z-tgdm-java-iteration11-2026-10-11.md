# P41z iteration-11 TGDM Java re-comparison — 2026-10-11

## Result

The selected Java plus TGDM sequence was replayed with the exact crop-ready
inputs saved from the P41z iteration-11 checkpoint. All six fields had zero
sign-changing voxels. Dice was 1.0 for each of the five extracted surfaces.
Maximum symmetric mean surface distance was `1.19e-6 mm`; maximum sampled
Hausdorff distance was `0.0109 mm`. All user-agreed surface gates passed
(Dice ≥0.97, mean distance ≤0.01 mm, Hausdorff ≤1 mm).

This supersedes the assumption that iteration-11 cortical support loss might
come from a Python TGDM mismatch. With the same saved source inputs, the Java
sequence reproduces the poor Python fields within tolerance. The saved case
still fails anatomy acceptance: only 0.1123% of label-100 GM is between the
surfaces and 59.51% of label-44 WM is inside the inner surface. On this crop,
the WM sign fraction falls from 99.38% in `initial_level_set` to 59.75% in
`initial_inner`, and only 0.208% of label-100 GM lies inside the final outer
surface. Dura GM itself is present at 66,330/67,677 label-100 voxels above
0.4, so the failure is not a complete loss of GM membership at native transfer.

The report is a held-fixed selected-method comparison. Its source signed-
distance seed is reused from the hash-checked Python checkpoint; tissue,
ACE/GVF, and TGDM inputs are serialized and replayed on the same 73×91×63 RAS
crop. It is not an independently generated Java workflow, a proof of correct
upstream inputs, or real-MRI validation. It shows that the selected TGDM
implementation itself matches Java on this poor-quality checkpoint.

## Reproduction

```sh
scripts/compare_plus_tgdm_saved_java.sh \
  build/p41z-case-normalization/iterations-0008-to-0011-surfaces/plus-surfaces \
  build/p41z-case-normalization/tgdm-java-iteration-0011 \
  docs/reports/p41z-tgdm-java-iteration11-2026-10-11.json \
  build/p41z-case-normalization/iterations-0008-to-0011-surfaces/synthetic-plus-checkpoint-surface-composition.json
```

The Java harness used **0.537 GB peak sampled RSS** under an 8 GB watchdog.
The complete per-field errors, input hashes, output hashes, and settings are
in the [machine-readable report](p41z-tgdm-java-iteration11-2026-10-11.json).
