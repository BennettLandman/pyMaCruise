# TOADS thinning-to-growth ordering state — 2026-10-05

## Source behavior and defect

The selected `MultipleToadSegmentation` constructor initializes its voxelwise
`ordering` array to zero (`MultipleToadSegmentation.java:314`). During
`propagateThinning`, each eligible popped voxel stores its queue priority at
that voxel (`:1962–1963`). The next `initGrowingLabels` adds the stored value
to each candidate's membership-ratio priority (`:493–494`). This field belongs
to the persistent segmentation object and therefore must survive thinning,
growth, and each subsequent loop pass.

The Python thinning result previously omitted this array, and the first
registered cycle created a fresh zero array for growth. That made the
individual growth kernel/source references pass while the composed Vandy
run diverged after the registration handoff.

## Implementation

`ToadsThinningPass`, `ToadsThinningQueueRun`, and individual thinning steps now
return the updated ordering array. The thinning APIs accept an initial array;
omitting it models the constructor's zero initialization. Eligible source
queue events write their popped priority into that array. Initial registered
growth now consumes the first thinning result; next-loop thinning starts from
the preceding growth result; later growth consumes the ordering returned by
that thinning. The continuing iteration-slice adapter follows the same order.

## Selected Java comparison

On the committed synthetic Vandy input (`36×45×34` native, `53×60×51`
working grid), the registered-loop output now compares with the pinned Java
reference as follows:

| State | Result |
|---|---:|
| First registered transform | Bit-identical |
| First grown segmentation | Exact, 0 of 162,180 voxels differ |
| Refreshed closest labels | Exact, 0 of 648,720 entries differ |
| Next thinning segmentation | Exact, 0 of 162,180 voxels differ |
| Iteration-2 centroids | Maximum absolute difference `7.1525574e-7` |
| Iteration-2 standard deviations | Maximum absolute difference `1.8114224e-7` |
| Iteration-2 memberships | MAE `2.6658711e-10`; maximum difference `4.3064356e-6`; zero values exceed `1e-5` |
| Iteration-2 objective | Same float32 value as Java (`43.742954`) |

These results match the selected synthetic Java run through this two-iteration
boundary; they do not establish convergence, the complete TOADS wrapper, the
full MaCRUISE workflow, or real-image validation.

## Focused validation

```bash
PYTHONPATH=src ./.venv/bin/pytest -q \
  tests/test_toads_vandy_iteration_boundary.py \
  tests/test_toads_registered_first_iteration_cycle.py \
  tests/test_toads_continuing_iteration_slice.py \
  tests/test_toads_default_11_continue_thinning.py \
  tests/test_toads_default_continue_thinning.py \
  tests/test_toads_next_default_membership.py \
  tests/test_toads_p18cyt_registered_iteration.py \
  tests/test_toads_p18cyv_next_iteration.py \
  tests/test_toads_p18cyw_exact_iteration.py \
  tests/test_toads_p18cyx_registered_cycle.py \
  tests/test_toads_p18cyz_next_loop_thinning.py \
  tests/test_toads_post_registration_next_thinning.py \
  tests/test_toads_thinning_queue.py \
  tests/test_toads_thinning_queue_step.py \
  tests/test_toads_thinning_queue_run.py
```

Result: **19 passed, 0 failed, 0 skipped in 34.27s**. This is focused
synthetic/source-reference validation. The full suite passed **1,202 tests, 0
failed, 0 skipped, and 11 expected RAS warnings in 1,187.29s** on implementation
commit `c8a33033e823d9a2fbce5170ee7ebd984008e801`, with peak sampled RSS
**13.539 GB** under the 24 GB watchdog. See the [full-suite report](full-pytest-2026-10-05-toads-ordering-carry.md).
