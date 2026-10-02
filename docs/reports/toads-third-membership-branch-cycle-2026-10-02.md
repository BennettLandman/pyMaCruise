# TOADS iteration-three membership through registration and wrapper branch (P18cybj)

## Ported behavior

Continue P18cybi's iteration-three exact membership result through the
selected source float32 energy ratio, rigid atlas registration, convergence
and iteration-limit checks, and terminal/continuing branch. The adapter
carries the immediately preceding thinning, grow ordering, volumes,
availability, relation distances, and transform into the next registration
call. A continuing branch refreshes its relations and shape priors and can
enter iteration-four thinning; a stopping branch transfers final labels to
the output classification map.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 26.91 seconds**.

The committed compact RAS fixture now runs through three membership updates and
three registration/control transitions. The test checks that registration
receives the exact third membership object, verifies energy ratio from the
selected float32 helper, and follows the selected branch. If that branch
continues, its next thinning arrays match a direct stage call. Registration
and control methods have selected Java-reference tests elsewhere; this is
synthetic composition, not combined Java wrapper execution or full TOADS
equivalence.

Full-suite evidence is pending.
