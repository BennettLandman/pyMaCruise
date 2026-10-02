# TOADS source-line accounting (started 2026-10-01; refreshed 2026-10-02)

## Current planning estimate

As re-reviewed on 2026-10-02 at implementation commit
`dcb707fec0425c3e28fe7b20bf2b2c5822a061be` with the full synthetic suite, the estimate remains
**3,400–4,800 source lines represented by Python
operations, out of 13,651 selected Java source lines (about 25–35%)**. Treat the
numerator as a coarse method-level estimate, not a measured line-coverage
result. The detailed acceptance rows and synthetic evidence define what is
represented; a bounded method or operation does not count its whole enclosing
class or make its wrapper complete.

The estimate covers the selected VandyCRUISE TOADS and inter-pass dura-removal
source below. It excludes visualization code and BET/SPECTRE/FreeSurfer because
the selected workflow does not port or require them. It also does not include
the separate LesionTOADS Java implementation; that work is tracked under P39.

## Fixed denominator

These are raw file line counts (`wc -l`) from the selected source checkout. The
denominator should be recomputed if the source revision or selected file set
changes.

| Selected source file | Lines |
|---|---:|
| `VandyCRUISE/src/edu/jhu/ece/iacl/plugins/classification/MedicAlgorithmToads08.java` | 383 |
| `VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java` | 1,156 |
| `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadAtlas.java` | 3,079 |
| `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java` | 3,324 |
| `VandyCRUISE/src/edu/jhmi/rad/medic/methods/InhomogeneityCorrection2.java` | 1,270 |
| `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadPreprocess.java` | 960 |
| `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java` | 3,184 |
| `VandyCRUISE/src/edu/jhu/ece/iacl/plugins/segmentation/skull_strip/MedicAlgorithmRemoveDura09.java` | 295 |
| **Total** | **13,651** |

## Numerator method and limitations

The numerator range is a planning estimate from the source methods and
operation families already represented by Python implementations and selected
synthetic Java/source comparisons (principally the P18 rows and their reports).
It counts selected algorithm operations, not equivalent Python LOC. One
concise NumPy/SciPy implementation may represent a longer Java loop, and a
partially ported class contributes only its represented operations. Rows that
pin fixed upstream inputs are evidence for the downstream operation only.

The percentage therefore helps track broad port scope, but it must not be
reported as executed Java-line coverage, complete method coverage, or workflow
equivalence. The broad range is intentional until a curated method-to-source-
line mapping is maintained. For each periodic plan update, refresh the range
from newly validated source operation groups and recompute the denominator
when the selected file set or source revisions change.

The TOADS startup composition through initial membership and first thinning
was extended through a synthetic selected 11-class atlas registration result
at implementation commit `dcb707fec0425c3e28fe7b20bf2b2c5822a061be`. This
latest adapter adds an actual registration-to-membership data handoff to the
validated compositions, but no additional Java algorithm method family; its
increment to the source-line numerator is therefore **0 lines**. It builds
on already-counted relation, deviation, image-variance, similarity,
shape-sampling, and membership operations. The source-line interval remains
3,400–4,800 of 13,651 lines (about 25–35%) pending a curated
method-to-source-line map. The full suite passed **1,060 tests**, with 11
expected RAS warnings in 1,862.37 seconds; see the [current full-suite
report](full-pytest-2026-10-02-registered-startup-membership.md).
