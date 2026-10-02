# TOADS iteration-four registration and terminal label transfer (P18cybl)

## Ported behavior

Continue the fourth exact membership result through selected rigid atlas
registration and the configured maximum-iteration stop. Verify the wrapper
returns no next thinning pass and transfers the current segmentation indices
to descriptor template labels while counting changed classification voxels.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 33.73 seconds**.

The compact synthetic RAS test checks membership identity into registration,
the terminal branch and absence of another thinning pass, then compares output
classification and changed-voxel count against `toads_classification_from_labels`.
Registration and iteration control have selected Java-reference component
evidence; this remains a composed synthetic run, not a combined Java wrapper
execution or full TOADS equivalence. Full-suite evidence is pending, so P18cybl
remains Started.
