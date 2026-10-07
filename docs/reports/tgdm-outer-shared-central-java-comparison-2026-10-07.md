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

The recorded comparison was rerun on code based on `ffb30e7`; after adding the
per-voxel iteration trace, the first difference appears in iteration 1 at crop
voxel `(73,76,73)`: Python clamps the crossing to `+0.1`, while Java accepts
`−0.01752247`. The direct Java topology harness and Python both reject the
initial neighborhood as a simple point. Earlier voxel updates in the same pass
may change that neighborhood before this point is visited; the current trace
does not yet capture the in-pass neighborhood or force at the exact update.
The focused TGDM tests passed **2/2 in 12.50s** before the trace hook was
expanded; the new per-iteration trace unit test passed **1/1 in 0.76s**.

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
simple-point checks agree, but the iteration trace finds a first-pass crossing
difference at `(73,76,73)`. Next, record the in-pass neighborhood, raw update
value, and topology result at that exact update, then compare the 33 final sign
disagreements. This remains
a synthetic source-comparison blocker; it does not establish workflow
equivalence or real-data behavior.
