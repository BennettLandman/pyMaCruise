# LesionTOADS registered centroid handoff — 2026-10-04

## Scope

P39zzzzzzbl carries the class centroids initialized on the normalized,
pre-registration crop through atlas alignment and the selected post-registration
recrop. It copies those values into the constructor's current/prior centroid
arrays. It also computes lesion and black-hole centroids from the recropped
normalized channels and applies the selected reference setters. It does not
compute the first post-constructor membership maps or a full LesionTOADS
iteration.

## Synthetic and selected-source evidence

The compact RAS T1/FLAIR intensity phantom has a checked-in Java
source-excerpt harness. The harness SHA-256 is
`f0c7d925ff96cc0414c04714fe0335ce03ff023f3a388b1e618a349d1262ae53`; its
generated `java_reference.json` SHA-256 is
`36cbef10107ca10d357efb994b2031522ee2725b7cb3b399290dada069fbc615`. Running
the documented `javac`/`java` commands reproduced the JSON byte-for-byte.
Python class-centroid bounds and values match the reference float32 values.

The composed RAS crop/recrop test verifies that the same pre-registration
centroid bits are copied to both constructor centroid arrays. It verifies that
the auxiliary centroids use the second crop, match direct source-derived
component calls, and are retained by reference by the Java-equivalent setters.
The test also checks wrapper stage order and input immutability. The
source-excerpt harness and existing component tests establish the centroid
formula and setter semantics; this is not a composite Java execution or a full
LesionTOADS output comparison.

Focused command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q \
  tests/test_lesion_toads_preprocessing.py \
  tests/test_lesion_loaded_atlas_initial_alignment.py \
  tests/test_lesion_registered_recrop.py \
  tests/test_lesion_registered_initial_thinning.py \
  tests/test_lesion_loaded_atlas_first_thinning.py \
  tests/test_lesion_initial_atlas_alignment.py \
  tests/test_lesion_initial_alignment_bounds.py
```

Result: **15 passed, 0 failed, 0 skipped in 14.11s** on implementation
commit `93fafc576ad98010a43ebe2e0afe33d74f5da756`.

## Full-suite evidence

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested commit: **`93fafc576ad98010a43ebe2e0afe33d74f5da756`**. Result:
**1,162 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 856.98s**.
Peak sampled RSS was **14.434 GB**; the watchdog did not terminate the run.

This is synthetic integration and selected-method source evidence. It does not
establish full LesionTOADS equivalence, MaCRUISEplus integration, or real-MRI
validation.
