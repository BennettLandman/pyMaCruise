# P18ap: TOADS inverse shape-prior energy

## Scope

This slice ports the positive-`priorCoefficient` atlas prior accumulation in
`MultipleToadSegmentation.computeExactMemberships`. For each interior voxel
and candidate class `k`, the selected source sums across `m != k` the product
of the other class's atlas shape-prior weight, the pairwise centroid
similarity coefficient, and the squared transformed shape value. It returns
both the addition to reciprocal-energy (`priorCoefficient * invprior`) and
the half-weighted objective-energy contribution
(`0.5 * priorCoefficient * invprior`). The source scans only interior
voxels; borders remain unchanged and are represented by zero contributions.

This is the atlas inverse-prior term only. It does not perform shape
registration or sampling, generate the transformed shapes, normalize
memberships, or validate a complete TOADS evolution.

## Synthetic Java formula reference

The deterministic 5×6×4×3 transformed-shape field is committed as a RAS
1 mm `.nii.gz` volume. Class shape-prior values and the similarity matrix
from the P18ao fixture are supplied to a Java formula harness with positive
coefficient `0.15`. It records every voxel/class inverse-energy and
objective-energy value. The selected source `MultipleToadSegmentation.java`
SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.

Regenerate the transformed-shape image with
`scripts/generate_toads_inverse_shape_prior_fixture.py` and the reference
output with
`scripts/generate_java_toads_inverse_shape_prior_reference.sh`.

## Validation

`tests/test_toads_inverse_shape_prior.py` compares both returned arrays
exactly as float32 with the Java reference. It verifies interior-only updates
and zero border contributions, no contribution for a zero coefficient,
exclusion of the candidate's own shape value, and invalid input dimensions.
The focused module passes **6 tests**.

The full configured suite passes **340 tests in 64.85 seconds** on
implementation commit
`1c70252edbb22adf544dac3df101d1e7df7fb98f`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This is formula
agreement for the synthetic inverse-prior contribution only.
