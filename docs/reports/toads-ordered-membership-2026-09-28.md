# TOADS ordered membership sweep

**Date:** 2026-09-28  
**Roadmap:** P18ak complete  
**Validated implementation/fixture commit:** `defd6c9c76b1c7a3f177eb2370e24d53621bd1b7`

## Selected source behavior

`src/pymacruise/toads_math.py::toads_ordered_membership_sweep` implements a
single interior membership sweep using the selected source loop order: x,
then y, then z, then class. For each candidate class it starts with a
precomputed data energy, computes the six-face squared-membership penalty
from the **current in-place volume**, adds the full smoothing term to
reciprocal energy and half the term to objective energy, applies the
direct/neighbor/prior relation factor and prior-only zero-shape rule, then
normalizes across classes. The objective accumulates in source voxel/class
order. Outer voxels are unchanged.

The selected `MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`; the
VandyCRUISE `Numerics.java` SHA-256 is
`13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb`.
`ToadsOrderedMembershipSweepReferenceHarness.java` is a standalone
source-formula excerpt for this selected sweep.

## Synthetic fixture and results

The comparison reuses the committed 7×6×5×4 RAS `.nii.gz` data-energy and
initial-membership volumes from `toads_direct_membership/`, plus shape priors
and closest-label inputs from `toads_shape_relations/`. The effective
smoothing coefficient is `0.017f`, smoothing priors are
`[0.5f, 0.8f, 1.1f, 0.25f]`, and `priorFactor` is `0.75f`.

All final memberships match the Java harness exactly after float32
conversion. The cumulative objective scalar also matches exactly. Five
focused tests pass, including a comparison proving the source in-place sweep
differs from a frozen-input synchronous smoothing update; input RAS geometry,
normalization, unchanged borders, and invalid inputs are checked.

## Scope limits

This is one small synthetic sweep over **precomputed** data energies. It does
not generate image likelihoods, PSF correction, atlas registration, shape
priors from the atlas, the local topology admissibility check, stopping or
iteration control, dura removal, or final TOADS segmentation equivalence.
The Java reference is a method-formula harness rather than a complete
Java/MIPAV pipeline output.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **317 passed in 51.09 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `defd6c9c76b1c7a3f177eb2370e24d53621bd1b7`
- Run date: 2026-09-28
