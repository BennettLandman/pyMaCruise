# Active LesionTOADS descriptor energy dispatch

This bounded adapter extends raw LesionTOADS data-energy construction from
the older four-class composition to the bundled 13-class descriptor. It uses
the descriptor's actual ordered topology labels and supports active OBJECT,
LESION, and MASK branches. The selected VandyCRUISE
`LesionToadSegmentation.computeMemberships` has its optimized-CSF PowerTable
implementation commented out; those descriptor entries therefore reach the
source fallback `num = outlier`. That behavior is represented explicitly.

The adapter does not implement membership smoothing, class/auxiliary energy
dispatch, inversion, relation gating, or complete loop integration. Its
scope is raw data terms. The committed Java formula harness pins the selected
source file and descriptor order and reproduces each class term bit-for-bit
for a deterministic interior voxel, including the active scalar outlier
fallback for optimized-CSF classes.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_active_class_energy.py tests/test_lesion_data_energy_maps.py
```

Result: **7 passed, 0 failed, 0 skipped in 0.55 seconds**. Regression checks
assert the bundled descriptor has 13 classes (9 OBJECT, 2 optimized CSF,
1 MASK, 1 LESION), compare float32 result bits for residual, mask, and fallback
terms, check untouched boundary faces and fielded OBJECT behavior, and reject
unknown types and missing required field maps.

## Full suite

Command:

```text
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result on tested commit `d60c99777b2a776b95a7898b4120918559c935b0`:
**1,146 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 804.55
seconds**. The watchdog completed without termination; peak sampled RSS was
15.764 GB under its 24 GB limit. This run used synthetic data only. It does
not establish the complete LesionTOADS iteration, source-equivalent workflow,
or real-MRI validation.
