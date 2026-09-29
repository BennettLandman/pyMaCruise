# Full pytest evidence: P39zzb three-auxiliary volume composition

- Date: 2026-09-29
- Implementation commit tested: `b44982bef1e9b87c34d597fa8c9fb9d008161908`
- Command: `.venv/bin/pytest -q`
- Result: **518 passed**, 2 expected RAS orientation warnings, in 87.97 seconds.
- Focused three-auxiliary and related sweep suite: 13 passed.

Two committed selected-Java synthetic references exercise LESION,
BLACKHOLE, and OUTMASK in both lesion-class precedence orderings. They compare
every interior class, lesion, black-hole, and outlier float32 membership bit,
including smoothing and prior terms, distance sentinels, zero/unit outlier
weights, shape gating, repeated auxiliary normalization, and boundary
preservation. This validates the selected combined formula, not the complete
iterative LesionTOADS workflow.
