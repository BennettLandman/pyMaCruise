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

## Corrected first-iteration state trace

On 2026-10-06, I repeated the selected Java trace after finding that an older
ignored diagnostic input writer serialized NumPy C order while MIPAV's
`ModelImage.importData` expects x-fastest storage. That exploratory C-order run
produced a different padded grid and is excluded from comparisons. The
corrected Java input was written in NumPy Fortran order and both implementations
then used the same `(53,60,51)` padded grid. The committed T1 values and geometry
were unchanged.

On that common grid, the Java and Python pre-update initialization membership
arrays are bit-identical (Java trace SHA-256
`b279944c954adc2240cd95768c79eba88b72259a8d337cd7c9742af3dc0afd0e`). The
startup assigned-class standard deviations, closest-label map, startup shape
samples, and first thinning segmentation also match exactly. The startup
membership update differs at 77,339 float entries, but the maximum absolute
difference is one float32 ULP (`1.1920929e-7`), the mean absolute difference is
`1.3466327e-10`, and no entry differs by more than `1e-6`.

For the first loop update, Python was seeded with the Java startup memberships,
Java iteration-1 centroids and thinned labels, Java closest-label state and
shape samples, and a unity correction field. The Python exact membership pass
then matches the Java iteration-1 output within one float32 ULP (maximum
`1.1920929e-7`, mean `8.5862491e-11`, zero values above `1e-6`); its objective
is the same float32 value (`50.3000908`). By contrast, the natural Python
iteration-1 membership array differs by at most `4.5299530e-6`, with 16 values
above `1e-6`. This confirms the selected first-update kernel on Java-seeded
state and points to amplification of the tiny natural-run roundoff difference
as the loop advances. It does not establish converged trajectory equivalence.

The Java traces were produced by `VandyCoreSyntheticRunExport` with a 6 GB JVM
heap cap, 12-iteration cap, and one registration update. The Python diagnostic
used `dump_python_toads_source_vandy_history.py` with iteration history enabled.
All raw trace files and temporary harness instrumentation remain under ignored
`build/java-reference/`; no production algorithm source changed. No pytest or
full-suite run was performed for this diagnostic, and the selected TOADS line
estimate is unchanged.

## Scope and next comparison

This is a state-seeded selected-method comparison on one compact synthetic
fixture. It is not whole-loop convergence equivalence, final MaCRUISE output
acceptance, a real-MRI test, or proof that the accumulated trajectory remains
within the requested segmentation tolerance on other inputs. With corrected
x-fastest Java input serialization, the earliest measured divergence remains
the natural iteration-1 exact update, while direct replay from Java state
matches within one ULP. The next TOADS comparison should trace source/Python
centroid calculation and the first carried state after that update, then
repeat the source-state replay at subsequent loop boundaries.

No full repository suite was run for this diagnostic; no source code changed.
