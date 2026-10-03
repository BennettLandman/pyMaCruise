# TOADS source-line accounting (started 2026-10-01; refreshed 2026-10-02)

## Current planning estimate

As re-reviewed on 2026-10-02 at full-suite commit
`f3b874b`, the estimate is
**about 3,464–4,864 source lines represented by Python
operations, out of 13,651 selected Java source lines (about 25–36%)**. Treat the
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

The latest validated additions represent approximately **30 additional
selected Java lines** for parsing and retrieving the fixed atlas Intensity
Variance Factor used by selected membership updates (`MultipleToadDeformableAtlas.java`
lines 270–279 and 1111–1131). They are separate from the previously counted
approximately **34 selected Java branch lines** in `MultipleToadSegmentation`:
the `useApprox`
initial thinning-state arm (source lines 516–537) and selected approximate
queue-event volume-update arms (source lines 1924–1928 and 1980–1985). These
are partial branches within methods already represented by the port; they do
not complete the approximate loop or establish its convergence. Those branch
lines are already included in the starting interval for this update. The
additional Intensity Variance Factor parser/accessor above advances that
interval by about 30, to 3,464–4,864. This remains a coarse estimate pending a
curated method-to-source-line map, not measured runtime code coverage.

The full suite passed **1,130 tests**, with 11 expected RAS warnings in
1,213.79 seconds on commit `f3b874b`; see the [kernel-optimization full-suite
report](full-pytest-2026-10-02-toads-kernel-optimization.md). This cycle adds
zero newly counted selected Java source lines: it batches already represented
data terms, vectorizes already represented relation factors, and adds a
Python-only state diagnostic. The suite also includes the restartable synthetic
base workflow, source-factor propagation, and a two-update Java reference
comparison. The base workflow still records both TOADS stages as capped and
nonconverged, and the 50-cycle synthetic loop misses its 1e-4 convergence
threshold. This is not whole-workflow equivalence or real-data validation. The
current approximate TOADS line estimate remains **3,464–4,864 of 13,651**.
