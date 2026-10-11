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
implementation (`a > 0 ? a : -a`). The previous objective is now recovered as
`3727.62158203125` from the hash-verified iteration-5 synthetic integration
summary (`f03d9f0d...`). It is the prior Python integration checkpoint value;
it has not been independently reproduced by a Java iteration-5 run.

For the 71×89×67 fixture, refreshed statistics, all membership/energy outputs,
the objective, and the ratio match Java bit-for-bit. Java/Python ratio bits are
both `1040630131`, objective bits are both `1162622407`, and both NaN counts
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
fixed membership pass. Registration and a selected wrapper transition are now
also compared below; growth on this state and iteration-to-iteration carry
remain open, and no full registered iteration or human-data equivalence is
claimed.

## Registration and nonterminal wrapper branch

The next selected boundary uses the exact iteration-6 membership update and
the 13 bundled RAS atlas priors. The P41z source-excerpt harness implements
the active single-rigid `LesionToadAtlas.registerAllShapes()` sequence using
the supplied `ParametricTransform` and `ImageFunctions`; Python and Java agree
bit-for-bit on each trial transform and energy, accepted/rejected trial counts,
energy change, lambda, final transform, and all 13 positive-support bounds and
counts. The initial transform and previous objective are taken from the
hash-verified iteration-5 integration checkpoint. That checkpoint is Python
integration evidence, not Java iteration-5 evidence.

The wrapper's `calcSegmentation()` energy-ratio/convergence branch is compared
to the selected Java expression: the ratio bits match, it continues, enters
growth, and refreshes relations. This test deliberately supplied an empty
skeleton so it isolates branch decisions; its empty growth queue is not growth
evidence.

During the Java comparison, the old vectorized affine-coordinate construction
produced a different positive-support count. Python bounds mapping now
evaluates image physical offsets, rotation, translation, template-resolution
division, and template-center addition in the float32 order used by
`ParametricTransform.imageToTemplate()`. The 71×89×67 comparison and a compact
nonidentity-transform phantom both match Java support bounds/counts exactly.

The held-fixed relation comparison on the saved thinned segmentation also
matches `computeRelations()` exactly for a relation map, relation-membership
table, and 14 diffusion passes. More importantly, the committed
`plus_iteration_6_thinning_result/` fixture supplies the true pre-growth
ordering and per-class volume after `initThinningLabels()` and
`propagateThinning()`. That result is byte-compared to the selected Java
thinning reference; its segmentation and skeleton match the P41z membership
fixture. These are valid inputs for the subsequent growth pass, not values
recovered from the later post-growth registration-loop outputs.

With the current iteration-6 membership update and that source-matched
pre-growth state, Python growth matches the selected Java harness byte for
byte. The comparison covers each popped heap item, topology decision, inserted
neighbor and priority, final segmentation, skeleton, availability, ordering,
and class volumes. `computeRelations()` after that grown segmentation also
matches the selected Java relation map and membership table. This validates
these held-fixed synthetic operations; it does not establish a complete
iteration workflow or human-MRI equivalence.

Focused evidence (no full suite):

- Registration + wrapper branch: **1 passed in 148.09s**, peak sampled RSS
  **0.499 GB** under a 4 GB watchdog.
- Growth + post-growth relation refresh: **1 passed in 116.56s**, peak sampled
  RSS **0.700 GB** under a 4 GB watchdog.
- Saved thinned-state relation comparison: **1 passed in 4.67s**.
- Compact transform-rounding/support-bound comparison: **1 passed in 1.29s**.

The long comparison is opt-in and reproducible with:

```sh
PYMACRUISE_RUN_P41Z_ITER6_REGISTRATION_BRANCH=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/pytest -q --tb=short \
  tests/test_p41z_iteration6_registration_branch_java.py::test_p41z_iteration6_registration_and_wrapper_branch_match_references
```

The growth and post-growth relation replay is opt-in and reproducible with:

```sh
PYMACRUISE_RUN_P41Z_ITER6_GROWTH=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/pytest -q --tb=short \
  tests/test_p41z_iteration6_registration_branch_java.py::test_p41z_iteration6_growth_and_postgrowth_relations_match_java
```

The relation and compact bounds comparisons are:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/pytest -q --tb=short \
  tests/test_p41z_iteration6_registration_branch_java.py::test_p41z_saved_thinned_state_relation_refresh_matches_selected_java \
  tests/test_p41z_iteration6_registration_branch_java.py::test_shape_bounds_transform_rounding_matches_java_on_compact_phantom
```
