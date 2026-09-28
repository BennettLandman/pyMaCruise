# P18ao: optimized-tissue centroid similarity

## Scope

This step ports the pairwise optimized-class similarity coefficient computed
in `MultipleToadSegmentation.computeExactMemberships`. It accumulates a
modality-weighted squared centroid distance normalized by the sum of each
class pair's standard deviations, then maps the distance through the selected
source formula with `priorScale=0.05f`. The matrix diagonal is explicitly
zeroed by the source.

This calculation produces class-pair coefficients only. It does not sample
the voxelwise shape prior, apply the coefficients to voxel energies, or
evolve memberships. The input classes are assumed to be optimized CSF, GM,
and WM classes.

## Synthetic Java formula reference

The 3×3 matrix fixture uses three modalities and classes. Its modality weights
come from the no-field output of the P18am Java reference, and its class
standard deviations come from the no-field P18an Java reference. The initial
centroids match the corresponding optimized-tissue synthetic case. The
reference records selected `MultipleToadSegmentation.java` SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.

Regenerate the formula output with
`scripts/generate_java_toads_optimized_similarity_reference.sh`.

## Validation

`tests/test_toads_optimized_similarity.py` checks the composed variance and
standard-deviation inputs against their upstream reference fixtures, then
compares all nine similarity coefficients exactly as float32 against the
Java formula harness. It checks a zero diagonal and pair symmetry for finite
positive standard deviations. A separate case preserves source IEEE behavior
when standard-deviation denominators are zero: different centroids yield
zero similarity and identical centroids yield NaN. The focused module passes
**2 tests**.

The full configured suite passes **334 tests in 51.06 seconds** on
implementation commit
`6f0a984c62f11c3b7bdfca57e154d7ff56215a9f`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This is formula-level
agreement for the synthetic class-pair calculation only.
