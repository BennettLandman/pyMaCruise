# Active LesionTOADS descriptor energy dispatch (started)

This bounded adapter begins extending LesionTOADS raw data-energy construction
from the older four-class composition to the bundled 13-class descriptor. It
uses the descriptor's actual ordered topology labels and supports active
OBJECT, LESION, and MASK branches. The selected VandyCRUISE
`LesionToadSegmentation.computeMemberships` has its optimized-CSF PowerTable
implementation commented out; those descriptor entries therefore reach the
source fallback `num = outlier`. That behavior is represented explicitly.

The adapter does not yet implement membership smoothing, class/auxiliary
energy dispatch, inversion, relation gating, or complete loop integration.
No Java formula harness is pinned for this slice yet, so its current evidence
is a synthetic Python contract check only.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_active_class_energy.py tests/test_lesion_data_energy_maps.py
```

Result: **6 passed, 0 failed, 0 skipped in 0.56 seconds**. The regression
checks assert the bundled descriptor has 13 classes (9 OBJECT, 2 optimized
CSF, 1 MASK, 1 LESION), verify residual/mask/fallback terms and untouched
boundaries, and reject unknown types and missing required field maps.

Next validation: add a harness pinned to the selected Java source and compare
float32 term results for all active descriptor families before marking the
plan row complete.
