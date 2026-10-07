# Isolated outer TGDM source comparison — 2026-10-07

## Scope and command

This synthetic comparison feeds the exact Java central level set from the
brain-resolution fixture, plus the same GM and ACE skeleton images, into the
selected Java and Python outer-surface calls. It isolates the outer update
from upstream central-field differences. The fixture is atlas-derived and
synthetic; no human MRI is included.

Regenerate the Java isolated output and run the Python comparison:

```bash
sh scripts/generate_java_tgdm_outer_stage_reference.sh
./.venv/bin/python scripts/compare_tgdm_outer_java_reference.py \
  --report docs/reports/tgdm-outer-shared-central-java-comparison-2026-10-07.json
```

The command exits nonzero because the Hausdorff gate is exceeded.

## Results

The pre-reinitialization fields have MAE `4.14e-6` and RMSE `9.10e-4`, but
contain 33 sign disagreements. The largest difference is at crop voxel
`(73,76,73)`: Python is `+0.1`, Java is `−0.779504`. The final output has Dice
`0.9999734`, symmetric mean surface distance `0.0003816 mm`, and Hausdorff
`1.257226 mm`; it misses the agreed `1 mm` Hausdorff limit. The JSON report
contains the voxelwise error counts and source profile details.

The selected Java subclasses pass distinct final distance caps to
`ProfileTGDM`, but that constructor does not assign the argument; the field
initializer remains 20. The Java run logs a 20 mm distance cap. The Python
profile retains this observed runtime behavior. A local experiment that used
the apparently intended 15 mm outer cap produced broad distance-field errors,
so that interpretation was discarded.

## Next work

The shared-central experiment localizes the remaining difference to outer
evolution before signed-distance reinitialization. The initial Python and Java
simple-point checks agree at the largest pre-reinitialization outlier, which
means later in-pass updates or update arithmetic change the neighborhood
before the crossing. Next, record the per-iteration Java/Python state around
the 33 differing sign voxels and match the first differing update. This remains
a synthetic source-comparison blocker; it does not establish workflow
equivalence or real-data behavior.
