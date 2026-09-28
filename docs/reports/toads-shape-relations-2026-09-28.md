# TOADS shape-relation factors

**Date:** 2026-09-28  
**Roadmap:** P18ai complete  
**Validated implementation/fixture commit:** `9e2e1d23ab2ecd325d382b871d4db98b7da89713`

## Selected source behavior

`src/pymacruise/toads_math.py::toads_shape_relation_factors` classifies each
candidate tissue from the four `closestLabel` slots used in
`MultipleToadSegmentation.computeExactMemberships`:

- a match in slots 0 or 1 is **direct**, factor `1`;
- otherwise a match in slots 2 or 3 is **neighbor**, factor
  `Numerics.bounded(2 * priorFactor * shape, 0.25, 1)`;
- otherwise it is **prior**, factor `min(1, 2 * priorFactor * shape)`;
- a prior-only candidate with shape exactly zero is marked suppressed, which
  corresponds to the source early branch that sets its membership to zero.

Direct classification has precedence over neighbor classification. The
selected `MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`; the
VandyCRUISE `Numerics.java` SHA-256 is
`13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb`. The
float `bounded` overload is `max(a, min(x, b))`.

## Synthetic fixture and results

The committed fixture contains 7×6×5×4 float32 shape priors and four-slot
integer closest labels, in RAS `.nii.gz` with a 1 mm isotropic affine. It
contains direct/neighbor collisions, empty (`-1`) slots, zero and fractional
shape values, and a prior factor chosen to exercise both clamp limits. Every
factor and prior-only zero-support decision matches the Java formula harness
exactly after float32 conversion. Eight focused tests pass, including manual
precedence cases and malformed inputs.

## Scope limits

This validates relation classification, scalar factors, and the zero-shape
early gate only. It does not apply factors to the reciprocal memberships,
compose normalization, determine the asynchronous update order, or validate
topology-constrained TOADS evolution. The fixture is a formula-level Java
reference, not a full segmentation comparison.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **305 passed in 59.40 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `9e2e1d23ab2ecd325d382b871d4db98b7da89713`
- Run date: 2026-09-28
