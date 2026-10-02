# LesionTOADS loaded-atlas first-thinning slice (2026-10-02)

## What now connects

Added `run_loaded_lesion_atlas_first_thinning()` in
`src/pymacruise/lesion_loaded_atlas_first_thinning.py`. It accepts the
descriptor-loaded, 13-class atlas alignment result and uses that exact initial
membership array as the membership input to the existing source-ordered
`initThinningLabels` and `propagateThinning` implementation. It checks the
class count and requires every supplied topology/frontier array to use the
loaded-membership grid.

The synthetic integration fixture extracts an 11³ patch from the committed
LesionTOADS shape atlas and creates normalized two-channel image values from
the corresponding descriptor means. It uses the resulting initial memberships
to produce a deterministic argmax segmentation and supplies explicit synthetic
skeleton, availability, ordering, and zero auxiliary lesion maps. The run
reaches and processes the actual thinning queue. No human imaging data is used.

## Source basis and limits

The wrapper call order is pinned in
`tests/data/synthetic/lesion_initial_atlas_alignment/java_reference.json` and
the source-pinned harness
`tests/java_reference/LesionInitialAtlasAlignmentReferenceHarness.java`.
That evidence covers `initialMemberships`, rigid registration setup, transform
handoff, and pyramid registration. The selected first-thinning order is
`initThinningLabels` followed by `propagateThinning`; those kernels have their
own Java reference fixtures in
`tests/test_lesion_toads_thinning_initialization.py` and
`tests/test_lesion_toads_thinning_pass.py`.

This new test validates wiring and successful execution on one atlas-derived
synthetic patch. It does not compare the combined 11³ result bit-for-bit with a
Java run. The synthetic argmax segmentation and explicit frontier/auxiliary
maps are integration inputs, not a port of their upstream construction. Thus
this is evidence that loaded atlas memberships can reach existing thinning
code, not evidence of an end-to-end LesionTOADS workflow or registration
accuracy.

## Validation

Implementation commit: `1896d7883c18085a0729cdfac6414d134ffbc531`.

Focused command:

```text
./.venv/bin/python -m pytest tests/test_lesion_loaded_atlas_first_thinning.py -q
```

Result: **2 passed in 2.36s**. The positive test reached a nonempty thinning
queue and produced events; the negative test confirmed that a grid mismatch
is rejected. `git diff --check` also passed before commit.

## Remaining work at this seam

Implement and validate source-derived segmentation/topology initialization,
distance maps, lesion/black-hole auxiliary memberships, and the transform of
shape priors used at this stage. Then compare the complete combined path to a
Java harness on a larger synthetic folded-cortex phantom before treating this
as more than a bounded integration slice.
