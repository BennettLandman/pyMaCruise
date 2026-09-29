# P18aq: ordered TOADS sweep with atlas prior energy

## Scope

This composition extends the source-ordered optimized-tissue membership
sweep with the atlas inverse-shape-prior terms from P18ap. At each interior
voxel, the sweep uses the current in-place membership field for six-face
smoothing, adds the atlas prior contribution to inverse-energy and its
half-weighted objective term, applies the direct/neighbor/prior relation
factor, converts energy to reciprocal membership, normalizes classes, then
accumulates the objective. Prior-only classes suppressed by zero shape skip
the prior calculation's effect, and volume borders remain unchanged.

This remains a small formula composition. It does not run topology
admissibility tests, generate the source data term, perform atlas sampling or
registration, iterate to convergence, or establish segmentation equivalence.

## Synthetic source-formula comparison

The Java harness reuses the committed 7×6×5×4 RAS 1 mm data-energy and
membership fixture and the transformed shape/closest-label fixture. It uses
four class shape-prior weights and a symmetric pairwise similarity matrix
computed from the selected source formula. The prior coefficient is `0.15`.
The harness records the final membership array and cumulative objective in
`tests/data/synthetic/toads_ordered_prior_sweep/java_reference.json`; its
source hash is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.

Regenerate the formula output with
`scripts/generate_java_toads_ordered_prior_sweep_reference.sh`.

## Validation

`tests/test_toads_ordered_prior_sweep.py` compares every normalized
membership and the source-ordered cumulative objective exactly as float32.
It also checks that a suppressed class ignores supplied prior terms, a zero
prior term reproduces the P18ak sweep, and image borders remain unchanged.
Together with the prior `tests/test_toads_ordered_membership.py` module, the
focused comparison/regression run passes **8 tests**.

The full configured suite passes **343 tests in 52.77 seconds** on
implementation commit
`0a53263aa8f70e8a5938ba53333f16f71af15299`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This is small-grid
formula-composition evidence, not a complete TOADS execution.
