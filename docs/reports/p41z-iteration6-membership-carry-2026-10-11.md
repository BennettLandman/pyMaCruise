# P41z iteration-6 membership carry checkpoint — 2026-10-11

## Source behavior added

On implementation commit `9e4e84a`, the ordered 13-class update preserves the
Java method-scope class-energy vector and lesion-energy scalar across its
x/y/z sweep. A skipped `shape[k] <= 0` class clears that class membership,
as `LesionToadSegmentation.computeMemberships()` does in its outer `else`
(source lines 1528–1530), while leaving the energy slot and lesion auxiliary
maps unchanged. The inner zero-shape branch around lines 1096–1104 is
unreachable because it is nested under `shape[k] > 0`.

The committed compact Java regression verifies an active lesion followed by
a skipped lesion class: class membership is cleared, while the lesion-energy
carry remains. The fixed P41z fixture has 222,045 interior class-10
nonpositive-shape voxels and no nonzero incoming class-10 membership among
them, so correcting the stale-membership behavior does not account for its
remaining energy discrepancy.

## Validation

- Focused voxel and compact Java regression: **8 passed** (commit `9e4e84a`).
- Full fixed-input P41z Java comparison: **not exact** on the same commit,
  under a 16 GB watchdog (peak sampled RSS 0.448 GB). The statistics refresh
  matches exactly; class memberships, class energies, lesion membership,
  black-hole membership, outlier membership, lesion relation, denominator,
  and black-hole energy all match bitwise. Java and Python each report zero
  NaNs.
- Remaining differences: lesion-energy values differ at 383,138 voxels
  (first `[1, 28, 29]`, max finite absolute error `0.72380799`); running
  objective differs at 388,088 voxels (first `[1, 33, 28]`, max finite
  absolute error `4.24145508`). Total objective float32 bits are Java
  `1162605034`, Python `1162622407`.

This is a held-fixed synthetic stage comparison, not full registered-loop
equivalence and not human-data validation. The exact remaining lesion-energy
source mismatch is unresolved; do not claim P41z membership parity. The
next diagnostic should inspect the first mismatch and its immediately prior
traversal voxel, including whether the lesion class is active and the raw,
smoothing, and prior contributions on both sides.
