# P41z iteration-6 Dura outputs: selected Java comparison

## Result

Using the committed, hash-provenanced RAS fixture in
[`tests/data/synthetic/plus_fine_grid_dura_iteration_6`](../../tests/data/synthetic/plus_fine_grid_dura_iteration_6/README.md),
the Java reference produced all four `LesionToadSegmentation.generateDuraRemovalOutputs`
channels exactly matching the saved Python outputs in float32 on all
423,373 voxels per channel:

| Output | Different voxels | Maximum absolute difference |
| --- | ---: | ---: |
| CSF | 0 | 0 |
| GM | 0 | 0 |
| WM fill | 0 | 0 |
| WM mask | 0 | 0 |

The arrays came from the state recorded as `completed_iteration: 6` (energy
ratio `0.13159732520580292`), not a newly generated approximation. The Java
reference input uses membership channels 0–12, lesion channel 13, black-hole
channel 14, and the matching segmentation. The registered checkpoint's
original SHA-256 values are retained in the fixture provenance. The selected
Java source SHA-256 is
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.

## Reproduction

```sh
.venv/bin/python -m pytest tests/test_plus_fine_grid_dura_java_reference.py -q
```

The focused test compiles the Java reference harness and requires a JDK. It
compares saved Python NIfTI outputs with the harness result; its machine report
is created in pytest's temporary directory. To create a standalone report:

```sh
.venv/bin/python scripts/compare_fine_grid_lesion_dura_java.py \
  --report build/p41z-dura-java-comparison.json
```

## Scope limits

This is exact component evidence for the Dura-output adapter on one synthetic
checkpoint. The harness mirrors the selected method's loop body, records the pinned
upstream source hash, and was separately audited against that source; it does
not instantiate the full Java segmentation class. The focused test uses only
committed NIfTI fixtures and the in-repository harness, so the original source
checkout is not required at test time. It does not establish upstream membership or segmentation
equivalence, later native-grid transfer equivalence, acceptable ribbon quality,
whole-workflow source equivalence, or real-human-MRI validation.
