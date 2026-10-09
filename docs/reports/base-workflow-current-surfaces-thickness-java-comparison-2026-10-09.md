# Corrected base workflow surface and thickness comparison against Java

## Result

The saved corrected seven-stage synthetic base run's `(172,216,164)` level
sets were passed to selected JIST `IsoSurfaceOnGrid.solveOriginal(...,18,0,
false)` and the production Python extractor. Each output mesh has identical
vertex and face counts, topology, and ordered triangle indices. Every vertex
maps bijectively, and the maximum corresponding-vertex displacement is
`7.1412e-6 mm` for each surface. This bounds continuous corresponding-triangle
distance below the agreed `0.01 mm` mean-position and `1 mm` Hausdorff limits.

| Surface | Vertices / faces | Ordered faces matched | Components / Euler | Max paired-vertex error |
|---|---:|---:|---:|---:|
| Inner | 115,518 / 231,032 | 231,032 / 231,032 | 1 / 2 | `7.1412e-6 mm` |
| Central | 120,726 / 241,448 | 241,448 / 241,448 | 1 / 2 | `7.1412e-6 mm` |
| Outer | 123,828 / 247,652 | 247,652 / 247,652 | 1 / 2 | `7.1412e-6 mm` |

The selected VandyCRUISE closest-point thickness wrapper was then run on the
same Java meshes and compared with PyMaCRUISE's production thickness API.
All seven physical fields matched exactly (maximum and mean absolute errors
were zero); central nearest-inner-vertex indices also matched exactly. The
Python API took `465.78 s` and peaked at `0.606 GB` process RSS. Java was
limited to a 2 GiB heap; the whole comparison ran under a 4 GB watchdog and
peaked at `0.534 GB` sampled RSS.

This is selected-source synthetic comparison on the current workflow's held-
fixed CRUISE level sets. It is not whole-workflow Java/MATLAB equivalence or
real-MRI validation. The Java mesh and thickness outputs and exact hashes are
retained under ignored `build/base-workflow-2026-10-09-retained-center/` paths.
The level-set scalar hashes, selected source hashes, topology, thickness
metrics, and runtime are in the [machine-readable surface
report](base-workflow-current-surfaces-thickness-java-comparison-2026-10-09.json).

## Reproduction

After running the corrected synthetic base workflow, run:

```sh
BASE=build/base-workflow-2026-10-09-retained-center
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 60 -- ./.venv/bin/python scripts/compare_workflow_jist_surfaces.py --core-dir "$BASE/stages/cruise-core/attempt-01" --work-dir "$BASE/java-jist-surfaces"
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 --interval-seconds 0.5 --report-interval-seconds 60 -- ./.venv/bin/python scripts/compare_jist_closest_point_thickness_brain_resolution.py --mesh-dir "$BASE/java-jist-surfaces/meshes" --levelset-dir "$BASE/stages/cruise-core/attempt-01" --work-dir "$BASE/java-jist-thickness" --report-json "$BASE/java-jist-thickness/comparison.json" --report-markdown "$BASE/java-jist-thickness/comparison.md" --title 'Corrected retained-center base workflow thickness source comparison'
```
