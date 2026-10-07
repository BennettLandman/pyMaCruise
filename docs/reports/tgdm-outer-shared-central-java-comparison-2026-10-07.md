# Brain-resolution TGDM Java source comparison — 2026-10-07

## Scope and commands

This comparison uses committed synthetic atlas-derived RAS `.nii.gz` inputs
on a `(166,208,156)` grid at 0.9375 mm spacing. No human MRI is included. It
compares the isolated outer solver with an identical Java central field and
then reruns the eight-stage TGDM core against the selected Java
`MaGenericTGDM`/`DistanceField` references.

```bash
sh scripts/generate_java_tgdm_outer_stage_reference.sh
./.venv/bin/python scripts/compare_tgdm_outer_java_reference.py \
  --report docs/reports/tgdm-outer-shared-central-java-comparison-2026-10-07.json
./.venv/bin/python scripts/compare_tgdm_core_java_reference.py \
  --report docs/reports/tgdm-core-brain-resolution-java-comparison-2026-10-07.json
./.venv/bin/pytest -q tests/test_tgdm.py
```

The harness records the source-ordered point value, proposed update, accepted
update, raw neighborhood and geodesic topology cube before Java's connectivity
checks mutate their input cube. Python records the same data for selected
points.

## Findings and correction

The first isolated update discrepancy was iteration 1 at crop voxel
`(73,76,73)`. Force values differed only by float precision, and both
implementations constructed the same topology cube. The Java connectivity
classes seed component searches from the center's 6- or 18-adjacent voxels.
They ignore components made only of unseeded corners. SciPy's full-cube label
count included such a corner component and rejected the update. The Python
simple-point check now counts only components reached from the source seed
neighbors, matching `SixConnected` and `EighteenConnected`.

The core comparison also exposed two tissue-sign voxels quantized one unit
below the Java result. Java evaluates `(1.0f - whiteMatterIso)` in float32;
Python had subtracted in double precision and cast afterward. The denominator
now follows the Java float32 operation order. A focused regression pins both
boundary values.

## Results

The isolated outer-surface comparison **passes**. Dice is `0.9999968`,
symmetric mean surface distance is `0.0000297 mm`, and Hausdorff distance is
`0.05685 mm`. The first traced update is accepted by both implementations;
the selected point neighborhood and topology cube match exactly. The
pre-reinitialization fields differ at four sign voxels, with MAE
`1.02e-6`; maximum absolute field difference is `0.05966`. After
reinitialization, the largest absolute field difference is `0.20266` while
surface distances remain inside the agreed limits. Full measurements and
point traces are in the linked JSON report.

The complete selected TGDM core comparison **passes** on the same crop. The
tissue-sign map is voxel-exact. Inner, central and outer masks have Dice
`1.0`, `1.0` and `0.9999887`; their symmetric mean surface distances are
`0.00000360`, `0.00000137` and `0.00014844 mm`; Hausdorff distances are
`0.00243`, `0.00329` and `0.19498 mm`. All are within the agreed thresholds
(Dice ≥0.97, mean surface distance ≤0.01 mm, Hausdorff ≤1 mm). The comparison
also reports all intermediate-field errors. The isolated outer comparison's
first update values differ by about `4e-10` from float arithmetic, but its
topology decision and accepted sign update match.

Focused TGDM module: **24 passed** (see the validated commit below). These are
synthetic selected-source comparisons and focused Python regression tests;
they do not establish full base-workflow equivalence or real-MRI behavior.

**Tested implementation commit:** recorded in the following documentation
update after committing the code and reference trace.
