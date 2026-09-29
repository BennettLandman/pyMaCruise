# P18at: local TOADS thinning transition

## Scope

This ports one already-popped interior queue item from
`MultipleToadSegmentation.propagateThinning`. It selects the first strict
maximum membership as the candidate tissue; examines six-face neighbors in
source order; prefers a directly neighboring candidate, otherwise chooses the
relation-connected neighbor with the highest transformed shape prior; calls
the previously validated multi-label 18/6 homeomorphic label-change gate; and
applies the strict `classVolume > minVolume` guard. For an eligible item it
removes the current class from the skeleton and decrements that class volume
in float32. It changes the segmentation label only if a proposed change passes
the topology gate. A popped item with no connected replacement still thins
the skeleton while leaving the label unchanged, matching source behavior.

The synthetic fixture contains nine neighborhood cases: unchanged best class,
exact volume threshold, stale skeleton item, no connected replacement, direct
accept/reject, relation-connected accept/reject, and a shape-prior tie resolved
by face iteration order. It records complete input arrays and expected state.
The harness uses the copied upstream `CriticalPointLUT` source and the bundled
`critical186LUT.raw.gz` reference asset.

Selected source SHA-256 values:

- `MultipleToadSegmentation.java`: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`
- `CriticalPointLUT.java`: `50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f`

Relevant selected source regions are `propagateThinning` (around lines
1863–2010), `bestMembership` (around line 1752), and
`isHomeomorphicSegmentation` (around lines 1554–1673). Regenerate the expected
states with `scripts/generate_java_toads_thinning_transition_reference.sh`.

## Validation and limits

`tests/test_toads_thinning_transition.py` compares every discrete output and
float32 class-volume value exactly for all nine cases. It also checks source
hashes, the strict threshold boundary, available flags, and accepted and
rejected topology outcomes. The focused thinning/topology/LUT modules pass
**15 tests**.

This validates one local state transition, not the priority queue, queue tie
behavior, whole-volume thinning initialization, repeated thinning/growing, the
two topology correction stages, atlas alignment, or final TOADS labels. Those
remain necessary before claiming TOADS segmentation equivalence.

The full configured suite passes **353 tests in 64.63 seconds** on
implementation commit `a04ae2ab66ece1ba1863522ef205e981bd8b5df2`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md).
