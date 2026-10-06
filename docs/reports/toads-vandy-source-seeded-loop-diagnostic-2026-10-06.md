# VandyCRUISE TOADS source-trace diagnosis — 2026-10-06

## Question

The compact synthetic TOADS loop diverged earlier than expected and stopped
one iteration later in the selected Java run. This investigation traced the
first membership differences to a source float32 coefficient-order mismatch,
corrected that mismatch, and repeated the stage and output comparisons.

## Inputs and source

- Synthetic input: the committed compact Vandy atlas T1 and brain-mask fixture
  under [`tests/data/synthetic/toads_downsampled_vandy_atlas/`](../../tests/data/synthetic/toads_downsampled_vandy_atlas/).
- Python driver and comparisons: [`test_toads_vandy_iteration_boundary.py`](../../tests/test_toads_vandy_iteration_boundary.py)
  and [`test_toads_effective_coefficients.py`](../../tests/test_toads_effective_coefficients.py).
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

## Cause and correction

On 2026-10-06, I repeated the selected Java trace after finding that an older
ignored diagnostic input writer serialized NumPy C order while MIPAV's
`ModelImage.importData` expects x-fastest storage. That exploratory C-order run
produced a different padded grid and is excluded from comparisons. The
corrected Java input was written in NumPy Fortran order and both implementations
then used the same `(53,60,51)` padded grid. The committed T1 values and geometry
were unchanged.

The selected Java constructor does not use the wrapper smoothing value
directly. `MultipleToadSegmentation` computes the effective coefficient as
`smooth / 6.0f / classes`, with each operation rounded to Java `float32`.
Python previously evaluated `0.2 / 6.0 / 11` in double precision and rounded
only once. For the 11-class profile these are adjacent float32 values:
Java-order `0x3b46980d`, prior Python value `0x3b46980c`. The prior coefficient
`2.0 / 11` already rounded to the same bit pattern in both implementations.

`toads_effective_smoothing` and `toads_effective_prior_coefficient` now encode
the source's operation order. The CLI and registered-loop reference drivers
use these helpers. On the corrected x-fastest common `(53,60,51)` grid, the
Java and Python pre-update memberships, startup statistics, closest-label map,
shape samples, and first thinning labels match exactly. With source-order
smoothing, the first three exact membership maps match every float32 bit;
iteration-2 centroids and standard deviations also match bit-for-bit. The
first observed natural-loop membership divergence is now iteration 4.

The corrected Python run reaches its configured energy threshold at iteration
11 (terminal ratio `7.6756049e-5`); the selected Java run reaches its 12-step
cap. The final class-label arrays are no longer byte-identical. Against the
committed Java final-label RAS NIfTI, all 11 per-label Dice scores are at least
`0.9975186` (the user's minimum is `0.97`), with four differing voxels. This is
compact synthetic acceptance evidence, not full workflow equivalence or real
MRI validation.

Focused validation: `./.venv/bin/pytest -q
tests/test_toads_effective_coefficients.py tests/test_toads_vandy_iteration_boundary.py`
passed **10 tests in 42.96s**. Base restart regression:
`./.venv/bin/pytest -q tests/test_base_workflow_restartable.py` passed **2 tests
in 37.02s**. The focused suite was rerun after cleanup and passed **10 tests in
41.32s**. Both validation commands exercised implementation commit
`9fe2ac6dd4c5adc173bda528b190f3be5df69b9f`. No full suite was run; the TOADS selected-source line estimate is
unchanged. Java reference-source files were not modified. Temporary raw trace
arrays and instrumentation remain under ignored `build/java-reference/`.

## Scope and next comparison

This is one selected VandyCRUISE synthetic run. It does not establish whole-loop
convergence equivalence or base/plus workflow acceptance. The next comparison
should trace the first differing iteration-4 state and later registration,
membership, and topology handoffs, then repeat the label Dice gate on a second
folded synthetic volume. Real-MRI acceptance remains deferred.
