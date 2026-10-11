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

- Focused voxel and compact Java regression: **8 passed** on implementation
  commit `9e4e84a` (with the Java oracle corrected in `733d19e`).
- Full fixed-input P41z Java comparison: **all outputs bitwise exact** on
  implementation commit `9e4e84a` with oracle/diagnostics commit `733d19e`.
  Statistics were exact first; memberships, class energies, lesion and
  black-hole auxiliaries, outlier, relation, denominator, lesion and
  black-hole energies, running objective, and total objective all match.
  Java and Python each report zero NaNs. The Java/Python total objective bits
  are both `1162622407`. The run used a 6 GB RSS watchdog and peaked at
  0.445 GB.
- The original apparent failure was in the Java reference harness: it
  initialized `eL=0` instead of source-faithfully assigning `eL=rawLesion`
  before smoothing and prior contributions. The full comparison now models
  the selected source assignment. At first mismatch coordinate `[1,28,29]`,
  a cropped traversal-prefix replay updated 811 preceding voxels; all six
  faces and raw, smoothing, prior, and final lesion-energy bits matched.
  The actual Vandy `PowerTable.java` lookup also matched.

This is a held-fixed synthetic stage comparison, not full registered-loop
equivalence and not human-data validation. The energy ratio/update and later
P41z registration/growth stages remain separately open.
