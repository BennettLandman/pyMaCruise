# TOADS shape-weighted membership update

**Date:** 2026-09-28  
**Roadmap:** P18aj complete  
**Validated implementation/fixture commit:** `7282405cdc3ead6ae00f495579240746f75aa7b4`

## Selected source behavior

`src/pymacruise/toads_math.py::toads_shape_weighted_membership_update`
composes the optimized-tissue direct relation update with the source shape
relation gate. On each interior voxel, a candidate class is first classified
as direct (factor 1), neighbor (bounded shape factor), or prior (shape factor).
A prior-only class with zero shape is set to zero before any reciprocal is
computed. All other classes use `1 / num` when `num > ZERO=1e-10f`, otherwise
the source `INF=1e10f`; each reciprocal is multiplied by its relation factor.
The resulting class values are accumulated in source order and normalized if
the denominator is positive. Outer memberships are untouched.

The selected `MultipleToadSegmentation.java` and `Numerics.java` source hashes
are documented with the underlying factor comparison in
[P18ai's report](toads-shape-relations-2026-09-28.md). The composition harness
retains the same source `closestLabel` ordering, float32 constants, branch
order, and formula-level reference limitation.

## Synthetic fixture and results

The comparison reuses committed synthetic RAS `.nii.gz` inputs from
`toads_direct_membership/` and `toads_shape_relations/`: data energy, initial
memberships, shape priors, and four-slot closest labels. This avoids duplicate
image fixture storage. The Java formula harness writes all updated
7×6×5×4 memberships to `tests/data/synthetic/toads_shape_weighted_membership/java_reference.json`.

Every membership matches the Java formula output exactly after float32
conversion. Tests cover all-unsupported candidates (denominator remains zero,
so updated candidates stay zero), mixed relation classes, normalization,
threshold energies, untouched borders, and invalid inputs. Seven focused tests
pass.

## Scope limits

This is a per-voxel composition over precomputed data energies with no
smoothing term. The complete source loop updates memberships asynchronously
and reads neighboring memberships while the sweep proceeds. That ordering,
topology admissibility, atlas registration, convergence, and full TOADS
segmentation remain open. This source-formula comparison is not a full Java
segmentation reference.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **312 passed in 50.92 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `7282405cdc3ead6ae00f495579240746f75aa7b4`
- Run date: 2026-09-28
