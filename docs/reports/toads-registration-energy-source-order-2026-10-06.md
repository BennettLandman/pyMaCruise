# TOADS registration energy operation-order correction — 2026-10-06

## Change

The selected VandyCRUISE `MultipleToadDeformableAtlas.registrationCost` computes each energy contribution as Java float32 operations in source order:

```text
membership * membership * prior * prior
```

The Python energy evaluators previously squared the membership and prior separately, then multiplied the squares. Those expressions are mathematically equivalent over real numbers but not bit-equivalent in float32. In the compact selected-source loop, the grouped version raised one candidate energy by one ULP. Python accepted that rigid update, while Java rejected it and kept the previous pose. That caused a pose difference at iteration 3 and changed later memberships.

The exact-membership, streaming, and vectorized energy paths now round after each left-associated multiplication. The selected-source registration pyramid harness had repeated the grouped expression; its Java test harness and committed numeric trace were corrected to follow `registrationCost` as well.

## Comparison fixture

- Evidence type: selected VandyCRUISE Java compact-synthetic reference and Python-to-Java comparison.
- Subject-like synthetic input: committed `(36,45,34)` T1 and matching mask; Vandy internally pads it to `(53,60,51)`.
- Registration source: `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`, SHA-256 `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`.
- Wrapper source: `AlgorithmMultipleToads.java`, SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`.
- Both runs used the committed 11-class Vandy atlas, one-level initialization, one registration update per iteration, source-order float32 smoothing, relation limit 3, and a 12-iteration maximum.

The Java output trace records the six rigid pose parameters after each registered loop cycle as exact uint32 float bit patterns in [`twelve_iteration_java_reference.json`](../../tests/data/synthetic/toads_downsampled_vandy_atlas/twelve_iteration_java_reference.json); the corresponding Java iteration memberships and final RAS NIfTI are committed alongside it.

## Results

- Registration transforms match Java bit-for-bit after cycles 1–11. At cycle 12 the largest parameter difference is `3.528594970703125e-05`.
- Membership arrays match bit-for-bit through iteration 10. At iteration 11, 279 of 1,783,980 membership values differ by more than `1e-6` (maximum absolute difference `8.64863395690918e-05`); at iteration 12, 3,717 exceed `1e-6` (maximum `0.0001436471939086914`). Membership argmax agrees at every voxel for all 12 iterations.
- Both source and Python runs reach the 12-iteration cap. Python's terminal energy ratio is `0.02333061955869198`; it is correctly reported as capped, not converged.
- Against the selected Java final RAS classification, Python has **0 changed voxels** and **Dice 1.0 for all 11 labels** on this compact fixture.

This demonstrates a source-faithful registration handoff and an exact final classification for this one compact synthetic case. It does not establish whole-workflow equivalence, natural convergence, folded-phantom coverage, or real-MRI validation.

## Validation

Focused command on implementation commit `2619a7698c40a3d43063a8b65dfa2f3adeb8c567`:

```bash
./.venv/bin/pytest -q \
  tests/test_toads_registration.py \
  tests/test_toads_registration_sampling.py \
  tests/test_toads_rigid_registration_update.py \
  tests/test_toads_vandy_iteration_boundary.py \
  tests/test_base_workflow_restartable.py
```

Result: **44 passed in 68.75s**, 0 failed, 0 skipped. The run covers float32 energy operation order, scalar/streaming/vectorized agreement, selected pyramid registration, 12 Vandy loop boundaries, final-label Dice, and base restartability. The source-reference harness was regenerated after correcting its multiplication order; compilation emitted one upstream `finalize()` deprecation warning. No full repository suite was run, so the selected-source line estimate at the top of `plan.md` is unchanged.
