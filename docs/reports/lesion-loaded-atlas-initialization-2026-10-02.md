# LesionTOADS loaded-atlas initialization slice (2026-10-02)

## What now connects

Added `initialize_lesion_atlas_alignment_from_descriptor()` in
`lesion_initial_atlas_alignment.py`. It loads the bundled 13-class atlas's
T1/FLAIR intensity means and RAS shape maps, checks that both use the same
class order, computes the selected initial memberships from caller-supplied
normalized T1_MPRAGE/FLAIR arrays, and forwards those memberships and all 13
shape priors into the existing selected rigid registration pyramid.

The synthetic integration test uses an 11³ crop of the actual committed atlas
to generate normalized T1/FLAIR values from the atlas descriptor. It enters
the selected membership and one coarse rigid-registration scale. The
registration optimizer is intentionally configured for zero updates: this
checks the real atlas loading and stage wiring without implying that a
preprocessing transform or optimizer has been validated on a human scan. The
membership volume from this composition is bit-identical to calling the
separately Java-referenced initial-membership API with those same loaded
centroids.

## Evidence

Focused validation:

```text
./.venv/bin/pytest -q tests/test_lesion_loaded_atlas_initial_alignment.py \
  tests/test_lesion_initial_atlas_alignment.py \
  tests/test_lesion_initial_memberships.py
4 passed in 3.32s
```

This uses the real descriptor and all 13 real shape maps. It does not add a new
Java harness for the descriptor-to-alignment wrapper; numerical behavior is
validated by the existing Java-referenced membership and rigid-pyramid stage
tests, and by exact identity to the already validated membership function.

## Current boundary

This advances through atlas-backed initial memberships and the selected
initial registration seam. It does not yet reach initial thinning or a
continuing iteration. Those require source-derived preprocessing and
transform setup plus matching topology, segmentation, distance, auxiliary
membership, and thinning-frontier state. Constructing any of those from the
atlas alone would invent algorithm behavior, so the new adapter keeps them as
explicit caller responsibilities. It also does not read or normalize a human
T1/FLAIR pair; normalization remains upstream.
