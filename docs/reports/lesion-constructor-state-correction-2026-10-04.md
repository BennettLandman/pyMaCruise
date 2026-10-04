# LesionTOADS constructor-state correction — 2026-10-04

## Finding and correction

The P39zzzzzzba Python adapter had initialized constructor `lesion_membership`
from a descriptor prior. That was not the selected VandyCRUISE constructor
behavior. `LesionToadSegmentation.java` initializes
`lesionMembership[x][y][z] = 0.0f` at line 470, while the initial image-based
membership values are produced separately by `LesionToadPreprocess`.

The previous Java formula harness also encoded a nonzero value, so its
reference did not verify that constructor assignment. The Python adapter now
always initializes this constructor buffer to float32 zero and no longer
accepts the misleading `lesion_prior` argument. The Java expression harness
and committed reference were corrected to the source's zero value. The
post-constructor descriptor prior remains available to the later membership
path; it is not copied into this state buffer.

## Validation

Focused command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_lesion_initial_segmentation_state.py tests/test_lesion_initial_thinning_composition.py tests/test_lesion_initial_loop_entry.py
```

Result: **6 passed, 0 failed, 0 skipped in 14.99s**. The selected Java source
file SHA-256 remains
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`; corrected
harness SHA-256 is
`68e93d57f3d31d54a9ea70e745bfaf1666d59b03b5f74c911ba88b77563a0249`.

Full repository command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result on implementation commit `44c2916`:
**1,159 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 862.74s**.
Peak sampled RSS was **14.140 GB**; the 24 GB watchdog completed normally.
This suite used synthetic tests and did not run the supplied human MRI.

This corrects one constructor-state field and its reference evidence. It does
not validate complete LesionTOADS, plus integration, or source-output
equivalence.
