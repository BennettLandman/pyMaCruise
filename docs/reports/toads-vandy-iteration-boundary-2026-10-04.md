# Vandy synthetic TOADS iteration-boundary comparison — 2026-10-04

## Scope

This compares the first two exact-membership updates and the next thinning
state in the selected VandyCRUISE TOADS implementation on the committed
`36×45×34` synthetic T1 and bundled atlas. It is a bounded Java/Python source
comparison, not full TOADS or MaCRUISE equivalence and not real-data
validation. No clinical image was used.

## Reproduction

The committed reference arrays are
`tests/data/synthetic/toads_downsampled_vandy_atlas/two_iteration_java_reference.npz`.
Their provenance and input/source hashes are in the neighboring
`two_iteration_java_reference.json`. The Java trace used `AlgorithmMultipleToads`
and `MultipleToadSegmentation` with one registration level, one registration
update per loop, padding 10, relation limit 3, and a two-iteration cap. Its
membership arrays use the padded `(53,60,51,11)` working grid. The Python
regression reconstructs the same setup with the production atlas loader:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_toads_vandy_iteration_boundary.py
```

Focused result after pinning the additional stage snapshots: **1 passed, 0
failed, 0 skipped in 10.99s**. The reference fixture stores x/y/z/class
float32 arrays and compact registration, growth, closest-label, and thinning
states in a compressed NPZ; tests do not require a Java runtime or source
checkout.

## Results

| State | Selected Java | Python | Comparison |
|---|---:|---:|---|
| Iteration 1 objective | `50.30009` | `50.30009079` | Difference below `2e-6` |
| Iteration 1 memberships | — | — | Mean absolute difference `5.65e-10`; maximum `3.78e-6`; 67 of 1,783,980 values exceed `1e-6` |
| Iteration 1 registered transform | `[0.00129835, 0.00129835, 0.00129835, -6.717452, -2.248702, -6.781790]` | `[0.00149813, 0.00149813, 0.00149813, -6.717252, -2.248502, -6.781591]` | All six parameters differ by approximately `0.0002`; this is the first observed split after membership update |
| Iteration 1 post-growth segmentation | — | — | 55 of 162,180 voxels differ |
| Iteration 1 closest-label state | — | — | 843 of 648,720 entries differ |
| Iteration 2 thinning segmentation | — | — | 42 of 162,180 voxels differ; all 11 class Dice scores are at least `0.99558` |
| Iteration 2 objective | `43.742954` | `43.74158859` | Absolute difference `0.0013654` |
| Iteration 2 memberships | — | — | Mean absolute difference `5.738e-6`; maximum `0.39403`; 10,403 of 1,783,980 values exceed `1e-5` |
| Iteration 2 energy ratio | `0.13944969` | `0.139480755` | The trajectories have split |

The first exact update is close to float32 rounding scale. The earliest
observed split is the next registration transform: all six parameters differ
from Java by about `0.0002`, followed by 55 post-growth label differences.
Closest-label state has 843 differing entries. This localizes the start of the
drift but does not identify its cause; registration initialization,
termination, and floating-point order remain candidates. The second
membership update has a large local maximum error even though its whole-tensor
mean remains small, and is not declared equivalent under a tolerance. The
intermediate thinning result meets the project’s 0.97 segmentation Dice floor,
but that does not validate later memberships, convergence, or final workflow
outputs.

## Source pins

- `AlgorithmMultipleToads.java` SHA-256:
  `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `MultipleToadSegmentation.java` SHA-256:
  `63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884`
- Compressed Java arrays SHA-256 (including the new iteration-1 state pins):
  `493ac4c3bfe0bcc0401a739b393eeb8474f29347494ec9f89b59bef6f90617be`

The reference uses the selected source run, while the Python result is a
production-loop synthetic run. Neither is a whole-workflow acceptance result.
