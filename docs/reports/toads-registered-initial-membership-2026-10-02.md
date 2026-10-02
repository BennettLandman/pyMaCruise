# P18cyat: registered startup to initial membership

## Scope

`compute_registered_toads_initial_membership` composes state already produced
by the selected atlas-startup adapter into the existing source-ordered
statistics and exact-membership pass. It reads the recropped RAS image,
initialized memberships and centroids, constructor segmentation, descriptor
class types, closest-label map, and transformed shape volumes from one
run-shaped object. It converts image axes from `(channel,x,y,z)` to
`(x,y,z,channel)` and leaves exponent, atlas/smoothing priors, tissue weights,
and effective coefficients explicit.

The integration fixture writes a synthetic default 11-class RAS atlas as
`.nii.gz` assets, runs `register_toads_from_atlas_assets`, then supplies its
actual `ToadsAtlasAlignmentRun` to the adapter. Separate negative tests use a
run-shaped state to reject missing data and a non-default descriptor class
order. This does not compare the combined composition to one Java wrapper
execution or cover correction-field fitting and later iterations.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_atlas_membership_integration.py tests/test_toads_registered_initial_membership.py tests/test_toads_initial_statistics_membership.py
```

Result: **6 passed, 0 failed, 0 skipped in 17.79 seconds**. Implementation
and test commit: `dcb707fec0425c3e28fe7b20bf2b2c5822a061be`.

Assertions compare the class deviations, image-variance weights, centroid
similarity, membership array, denominator, and objective against direct calls
using the corresponding synthetic upstream arrays. Negative cases reject a
missing run state and a non-default descriptor class order.

## Full-suite evidence

The full command `./.venv/bin/pytest -q` passed **1,060 tests, 0 failed,
0 skipped, and 11 expected RAS-orientation warnings in 1,862.37 seconds
(31:02.37)** on the same commit. See the [full-suite report](full-pytest-2026-10-02-registered-startup-membership.md).

This cycle adds validated synthetic integration coverage, but no separately
countable Java algorithm method family. The selected TOADS source-line
estimate therefore records **0 newly countable source lines** for this cycle;
the adapter connects operations already represented in the estimate.
