# Full pytest after LesionTOADS descriptor statistics/membership composition

Command:

```text
./.venv/bin/pytest -q
```

Tested implementation commit: `936ba20da8ea6ab3ea3bdb238cc1ac100790c0dc`.

Result: **1,153 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
828.21 seconds (13:48)**. Peak RSS was not sampled during this run.

The suite includes the optional FLAIR-weighted descriptor lesion-energy path,
invalid channel/exponent checks, the 5×5×5 synthetic statistics-to-membership
composition, and repository-wide regressions. No human MRI data was run.
Existing Java-matched statistics, energy, and membership components are
composed here; this does not add a new Java reference for the complete adapter.
It validates a fixed-state transition, not iterative LesionTOADS, complete
MaCRUISEplus, workflow equivalence, or real-data behavior. No selected base
TOADS source lines were newly counted; the coarse estimate remains
3,464–4,864 of 13,651 (25–36%).
