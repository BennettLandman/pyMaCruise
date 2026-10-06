# VandyCRUISE source-seeded TOADS iteration 4 — 2026-10-06

## Question

The compact 12-iteration selected-Java comparison has matching terminal labels
but diverging intermediate memberships and energy ratios. This diagnostic
replays iteration 4 from the same source state to determine whether the first
degree-one correction fit or the exact membership update introduces the
observed drift.

## Inputs and source

- Synthetic input: the committed compact Vandy atlas T1 and brain-mask fixture
  under [`tests/data/synthetic/toads_downsampled_vandy_atlas/`](../../tests/data/synthetic/toads_downsampled_vandy_atlas/).
- Python driver and comparison: `_registered_two_iteration_run(maximum_iterations=4)` in
  [`test_toads_vandy_iteration_boundary.py`](../../tests/test_toads_vandy_iteration_boundary.py).
- Selected source: VandyCRUISE `AlgorithmMultipleToads.java`, SHA-256
  `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`, and
  `MultipleToadSegmentation.java`, SHA-256
  `63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- Run settings: the subject input is the committed `(36,45,34)`
  fixture; Java internally crops/pads it to `(53,60,51)`. The selected atlas
  has 11 classes, uses one registration level and one registration update,
  padding 10, relation limit 3, and an iteration cap of 4 for this diagnostic.
- The existing committed Java arrays in `twelve_iteration_java_memberships.npz`
  provide the iteration-3 and iteration-4 exact-membership states. Temporary
  Java diagnostic state dumps were written only under ignored `build/`.

## Source-seeded replay findings

The natural Python and Java runs have a small incoming state difference by
iteration 3: the largest membership difference is `2.9206276e-6` (96 values
exceed `1e-6`). The Python/Java label map after iteration-4 thinning and the
four-slot closest-label relation map are identical. Their transformed shape
priors for the next exact pass are bit-identical.

When Python's degree-one refit and next exact pass start from Java's
iteration-4 input memberships, centroids, and field state:

- the fitted degree-one field is bit-identical to Java;
- refreshed centroids, assigned-class standard deviations, and normalized
  image-variance weights are bit-identical to Java;
- the exact membership output differs from Java by at most
  `1.1920929e-7` (one float32 ULP), with mean absolute difference
  `8.8864381e-11`.

Using either the source-dumped Java field or the Python-refitted field from
that same state gives the same reported membership distance. The replay
therefore localizes the later natural-run drift to already-carried state
differences; it does not find a new mismatch in the iteration-4 field or
membership kernel.

## Scope and next comparison

This is a state-seeded selected-method comparison on one compact synthetic
fixture. It is not whole-loop convergence equivalence, final MaCRUISE output
acceptance, a real-MRI test, or proof that the accumulated trajectory remains
within the requested segmentation tolerance on other inputs. The earliest
measurable natural-run divergence appears at the first exact update: startup
memberships differ by at most one ULP, while iteration-1 memberships differ
by at most `4.5299530e-6` at 16 values. The next TOADS comparison should seed
the first loop update from identical Java state and compare its thinning
priority/state and exact-membership stages separately.

No full repository suite was run for this diagnostic; no source code changed.
