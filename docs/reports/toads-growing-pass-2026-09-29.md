# P18bc — initialized TOADS growing pass

## Scope

`toads_propagate_growing` composes selected source `initGrowingLabels` queue
initialization with `propagateGrowing` events until the Java-compatible
min-heap is empty. Each returned event carries the popped item, topology
outcome, and any newly scheduled face-neighbors. Final output retains the
segmentation, skeleton, availability, ordering, class volumes, and residual
queue state.

## Synthetic Java-reference comparison

The 9×9×9 four-class phantom contains a folded stair-stepped tissue band, a
compact second tissue region, and an edge-adjacent label region. Two skeleton
seeds start competing fronts; memberships, ordering, and initial class volumes
are deterministic float32 inputs. The selected actual 18/6 LUT and Java
`BinaryHeap4D` are used by the Java harness.

All **448 queue events** match in exact event order and float32 priority bits.
The trace has 337 accepted topology events, 60 rejected events, and 51 stale
events. Every event's topology outcome and inserted-neighbor list match, as do
the final segmentation, skeleton, availability, ordering bits, volume bits,
and empty residual queue. Focused tests for this pass and its component
initialization/event stages pass **3 tests**. The reference generator is
`scripts/generate_java_toads_growing_pass_reference.sh`; the formula and
phantom construction are documented in
`tests/data/synthetic/toads_growing_pass/README.md`.

## Evidence and limits

The complete configured pytest suite passed **364 tests in 74.98 seconds** on
implementation commit `90bb1424ece7378e206748ce7555489474cd34b7`; see
[`full-pytest-2026-09-29.md`](full-pytest-2026-09-29.md). This source-slice
pass validates queue exhaustion and local topology-gated growth on synthetic
labels only. Atlas registration, the rest of TOADS iterations, dura removal,
full workflow equivalence, and actual imaging data remain outside this step.
