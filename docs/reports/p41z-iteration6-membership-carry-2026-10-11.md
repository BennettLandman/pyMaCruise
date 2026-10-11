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
equivalence and not human-data validation.

## Wrapper energy-ratio boundary

The comparator now continues through the wrapper's next source expression in
`AlgorithmLesionToads.calcSegmentation()` (membership objective update and
`Eratio = Numerics.abs(2.0f*(Eprev-energy)/(Eprev+energy))`, around source
lines 687–706). `Numerics.abs(float)` is the source comparison-and-negation
implementation (`a > 0 ? a : -a`). The fixture does not include the true
previous objective, so both Java and Python receive an explicit synthetic
`Eprev=3.25`; this isolates and validates the loop expression and its
float32 wiring without claiming an iteration-5-to-6 objective history.

For the 71×89×67 fixture, refreshed statistics, all membership/energy outputs,
the objective, and the ratio match Java bit-for-bit. Java/Python ratio bits are
both `1073708481`, objective bits are both `1162622407`, and both NaN counts
are zero. The fixed-volume comparator completed under a 4 GB watchdog with
0.445 GB peak sampled RSS. The compact regression is:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/pytest -q --tb=short \
  tests/test_p41z_iteration6_membership_java.py::test_java_membership_energy_carries_from_active_lesion_to_skipped_support
```

Result: **1 passed**. The fixed-volume comparison is reproducible with:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python scripts/compare_p41z_iteration6_membership_java.py
```

This closes only the source-faithful energy-ratio expression after the held-
fixed membership pass. P41z registered atlas update and subsequent growth,
relation refresh, and iteration-to-iteration state carry remain open; no
full registered iteration or human-data equivalence is claimed.
